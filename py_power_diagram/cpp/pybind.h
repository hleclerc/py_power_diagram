// #include <pybind11/functional.h>
// #include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>

#include "../../ext/power_diagram/src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"

#include "../../ext/power_diagram/src/PowerDiagram/Visitors/ZGrid.h"

#include "../../ext/power_diagram/src/PowerDiagram/get_integrations.h"
#include "../../ext/power_diagram/src/PowerDiagram/VtkOutput.h"

namespace py = pybind11;

struct PyPc {
    static constexpr int nb_bits_per_axis = 31;
    static constexpr int allow_ball_cut   = 0;
    static constexpr int dim              = PD_DIM;
    using                TI               = std::size_t;
    using                TF               = PD_TYPE; // boost::multiprecision::mpfr_float_100
};

struct PyZGrid {
    using Grid = PowerDiagram::Visitor::ZGrid<PyPc>;
    using Pt   = typename Grid::Pt;
    using TF   = typename Grid::TF;

    PyZGrid( int max_dirac_per_cell ) : grid( max_dirac_per_cell ) {
    }

    void update( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, bool positions_have_changed, bool weights_have_changed ) {
        auto buf_positions = positions.request();
        auto buf_weights = weights.request();
        if ( buf_positions.shape[ 1 ] != PyPc::dim )
            throw pybind11::value_error( "dim does not correspond to shape[ 1 ] of positions" );

        grid.update(
            reinterpret_cast<const Pt *>( buf_positions.ptr ),
            reinterpret_cast<const TF *>( buf_weights.ptr ),
            positions.shape( 0 ),
            positions_have_changed,
            weights_have_changed
        );
    }

    void display_vtk( const char *filename ) {
        VtkOutput<0,TF> vo;
        grid.display( vo );
        vo.save( filename );
    }

    Grid grid;
};

struct PyConvexPolyhedraAssembly {
    using TB = PowerDiagram::Bounds::ConvexPolyhedronAssembly<PyPc>;
    using TF = PD_TYPE; // boost::multiprecision::mpfr_float_100
    using Pt = TB::Pt;

    PyConvexPolyhedraAssembly() {
    }

    void add_box( py::array_t<PD_TYPE> &min_pos, py::array_t<PD_TYPE> &max_pos, PD_TYPE coeff, std::size_t cut_id ) {
        auto buf_min_pos = min_pos.request(); auto ptr_min_pos = (PD_TYPE *)buf_min_pos.ptr;
        auto buf_max_pos = max_pos.request(); auto ptr_max_pos = (PD_TYPE *)buf_max_pos.ptr;
        if ( min_pos.size() != PyPc::dim )
            throw pybind11::value_error( "wrong dimensions for min_pos" );
        if ( max_pos.size() != PyPc::dim )
            throw pybind11::value_error( "wrong dimensions for max_pos" );
        bounds.add_box( ptr_min_pos, ptr_max_pos, coeff, cut_id );
    }

    void display_boundaries_vtk( const char *filename ) {
        VtkOutput<1,TF> vo;
        bounds.display_boundaries( vo );
        vo.save( filename );
    }

    TB bounds;
};

py::array_t<PD_TYPE> integration( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid, const std::string &func ) {
    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    py::array_t<PD_TYPE> res;
    res.resize( { positions.shape( 0 ) } );
    auto buf_res = res.request();
    auto ptr_res = (PD_TYPE *)buf_res.ptr;

    if ( func == "1" || func == "unit" )
        PowerDiagram::integration( ptr_res, py_grid.grid, domain.bounds, ptr_positions, ptr_weights, positions.shape( 0 ), FunctionEnum::Unit    () );
    else if ( func == "exp(-r**2)" || func == "exp(-r^2)" )
        PowerDiagram::integration( ptr_res, py_grid.grid, domain.bounds, ptr_positions, ptr_weights, positions.shape( 0 ), FunctionEnum::Gaussian() );
    else if ( func == "r**2" || func == "r^2" )
        PowerDiagram::integration( ptr_res, py_grid.grid, domain.bounds, ptr_positions, ptr_weights, positions.shape( 0 ), FunctionEnum::R2      () );
    else
        throw pybind11::value_error( "unknown function type" );

    return res;
}

void display_vtk( const char *filename, py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid ) {
    VtkOutput<1> vtk_output( { "weight" } );

    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    py_grid.grid.for_each_laguerre_cell(
        [&]( auto &lc, std::size_t num_dirac_0 ) {
            domain.bounds.for_each_intersection( lc, [&]( auto &cp, auto space_func ) {
                cp.display( vtk_output, { ptr_weights[ num_dirac_0 ] } );
            } );
        }, domain.bounds.englobing_convex_polyhedron(),
        ptr_positions,
        ptr_weights,
        positions.shape( 0 )
    );

    vtk_output.save( filename );
}

PYBIND11_MODULE( PD_MODULE_NAME, m ) {
    m.doc() = "Power diagram tools";

    py::class_<PyZGrid>( m, "ZGrid" )
        .def( py::init<int>()                                                             , "" )
        .def( "update"                , &PyZGrid::update                                  , "" )
        .def( "display_vtk"           , &PyZGrid::display_vtk                             , "" )
    ;

    py::class_<PyConvexPolyhedraAssembly>( m, "ConvexPolyhedraAssembly" )
        .def( py::init<>()                                                                , "" )
        .def( "add_box"               , &PyConvexPolyhedraAssembly::add_box               , "" )
        .def( "display_boundaries_vtk", &PyConvexPolyhedraAssembly::display_boundaries_vtk, "" )
    ;    

    m.def( "integration", &integration );
    m.def( "display_vtk", &display_vtk );
}

