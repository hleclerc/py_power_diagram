// #include <pybind11/functional.h>
// #include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>

#include "../ext/power_diagram/src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../ext/power_diagram/src/PowerDiagram/Visitors/ZGrid.h"
#include "../ext/power_diagram/src/PowerDiagram/VtkOutput.h"

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

    PyZGrid( int max_dirac_per_cell ) /*: grid( max_dirac_per_cell ) */{
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

    void display_vtk( const char *vtk ) {
        VtkOutput<0,TF> vo;
        grid.display( vo );
        vo.save( vtk );
    }

    Grid grid;
};

py::array_t<PD_TYPE> get_measures( const py::array_t<PD_TYPE> &positions, const py::array_t<PD_TYPE> &weights ) {
    std::cout << positions.shape( 0 ) << std::endl;
    std::cout << positions.shape( 1 ) << std::endl;

    py::array_t<PD_TYPE> res;
    res.resize( { 10ul, 3ul } );
    //        //        auto buf_res = res.request();
    //        //        auto ptr_res = (PD_TYPE *)buf_res.ptr;
    //        //        for( size_t i = 0, o = 0; i < pb.size(); ++i ) {
    //        //            ptr_res[ o++ ] = pb[ i ].x;
    //        //            ptr_res[ o++ ] = pb[ i ].y;
    //        //        }
    return res;
}

void display_vtk( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, const char *filename, py::object &py_grid ) {
    PyZGrid g = py_grid.cast<PyZGrid>();
        // VtkOutput<1> vtk_output( { "weight" } );

        // auto buf_positions = positions.request();
        // auto buf_weights = weights.request();

        // auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
        // auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

        // using Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<PyPc>;
        // Bounds bounds;
        // bounds.add_box( { 0, 0 }, { 1, 1 }, 1.0, -1 );

        // py_grid.grid.for_each_laguerre_cell(
        //     [&]( auto &lc, std::size_t num_dirac_0 ) {
        //         lc.display( vtk_output, { ptr_weights[ num_dirac_0 ] } );
        //     }, bounds.englobing_convex_polyhedron(),
        //     ptr_positions,
        //     ptr_weights,
        //     positions.shape( 0 )
        // );

        // vtk_output.save( filename );
}

PYBIND11_MODULE( PD_MODULE_NAME, m ) {
    m.doc() = "Power diagram tools";

    py::class_<PyZGrid>( m, "ZGrid" )
        .def( py::init<int>()                     , "" )
        .def( "update", &PyZGrid::update          , "" )
        .def( "display_vtk", &PyZGrid::display_vtk, "" )
    ;

    m.def( "get_measures", &get_measures );
    m.def( "display_vtk" , &display_vtk );
}

