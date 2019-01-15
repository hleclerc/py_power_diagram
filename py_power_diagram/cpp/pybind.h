// #include <pybind11/functional.h>
// #include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>

#include "../../ext/power_diagram/src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"

#include "../../ext/power_diagram/src/PowerDiagram/Visitors/ZGrid.h"

#include "../../ext/power_diagram/src/PowerDiagram/get_der_integrals_wrt_weights.h"
#include "../../ext/power_diagram/src/PowerDiagram/get_integrals.h"
#include "../../ext/power_diagram/src/PowerDiagram/get_centroids.h"
#include "../../ext/power_diagram/src/PowerDiagram/VtkOutput.h"

namespace py = pybind11;

template<class FU>
void find_radial_func( const std::string &func, const FU &fu ) {
    if ( func == "1" || func == "unit" ) {
        fu( FunctionEnum::Unit() );
        return;
    }

    if ( func == "exp(-r**2)" || func == "exp(-r^2)" ) {
        fu( FunctionEnum::Gaussian() );
        return;
    }

    if ( func == "r**2" || func == "r^2" ) {
        fu( FunctionEnum::R2() );
        return;
    }

    if ( func == "in_ball(weight**0.5)" ) {
        fu( FunctionEnum::InBallW05() );
        return;
    }

    throw pybind11::value_error( "unknown function type" );
}

struct PyPc {
    static constexpr int nb_bits_per_axis = 31;
    static constexpr int allow_ball_cut   = 1;
    static constexpr int dim              = PD_DIM;
    using                TI               = std::size_t;
    using                TF               = PD_TYPE; // boost::multiprecision::mpfr_float_100
};

struct PyZGrid {
    using Grid = PowerDiagram::Visitor::ZGrid<PyPc>;
    using Pt   = typename Grid::Pt;
    using TF   = typename Grid::TF;

    PyZGrid( int max_dirac_per_cell, PD_TYPE max_delta_weight_per_grid ) : grid( max_dirac_per_cell, max_delta_weight_per_grid ) {
    }

    void update( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, bool positions_have_changed, bool weights_have_changed, std::string radial_func ) {
        auto buf_positions = positions.request();
        auto buf_weights = weights.request();
        if ( buf_positions.shape[ 1 ] != PyPc::dim )
            throw pybind11::value_error( "dim does not correspond to shape[ 1 ] of positions" );

        grid.update(
            reinterpret_cast<const Pt *>( buf_positions.ptr ),
            reinterpret_cast<const TF *>( buf_weights.ptr ),
            positions.shape( 0 ),
            positions_have_changed,
            weights_have_changed,
            radial_func == "in_ball(weight**0.5)"
        );
    }

    void display_vtk( const char *filename ) {
        VtkOutput<1,TF> vo({ "num" });
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

    void add_convex_polyhedron( py::array_t<PD_TYPE> &positions_and_normals, PD_TYPE coeff, std::size_t cut_id ) {
        auto buf_pan = positions_and_normals.request(); auto ptr_pan = (PD_TYPE *)buf_pan.ptr;
        if ( positions_and_normals.shape( 1 ) != 2 * PyPc::dim )
            throw pybind11::value_error( "wrong dimensions for positions_and_normals" );
        std::vector<Pt> positions, normals;
        for( pybind11::ssize_t i = 0; i < positions_and_normals.shape( 0 ); ++i ) {
            positions.push_back( ptr_pan + PyPc::dim * ( 2 * i + 0 ) );
            normals  .push_back( ptr_pan + PyPc::dim * ( 2 * i + 1 ) );
        }
        bounds.add_convex_polyhedron( positions, normals, coeff, cut_id );
    }

    void display_boundaries_vtk( const char *filename ) {
        VtkOutput<1,TF> vo;
        bounds.display_boundaries( vo );
        vo.save( filename );
    }

    PD_TYPE coeff_at( py::array_t<PD_TYPE> &point ) {
        auto buf_point = point.request(); auto ptr_buf_point = (PD_TYPE *)buf_point.ptr;
        if ( point.size() != PyPc::dim )
            throw pybind11::value_error( "wrong dimensions for point" );
        return bounds.coeff_at( ptr_buf_point );
    }

    TB bounds;
};

py::array_t<PD_TYPE> get_integrals( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid, const std::string &func ) {
    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    py::array_t<PD_TYPE> res;
    res.resize( { positions.shape( 0 ) } );
    auto buf_res = res.request();
    auto ptr_res = (PD_TYPE *)buf_res.ptr;

    find_radial_func( func, [&]( auto ft ) {
        PowerDiagram::get_integrals( ptr_res, py_grid.grid, domain.bounds, ptr_positions, ptr_weights, positions.shape( 0 ), ft );
    } );

    return res;
}

py::array_t<PD_TYPE> get_centroids( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid, const std::string &func ) {
    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    py::array_t<PD_TYPE> res;
    res.resize( { positions.shape( 0 ), (py::ssize_t)PyPc::dim } );
    auto buf_res = res.request();
    auto ptr_res = (PD_TYPE *)buf_res.ptr;

    find_radial_func( func, [&]( auto ft ) {
        PowerDiagram::get_centroids( py_grid.grid, domain.bounds, ptr_positions, ptr_weights, positions.shape( 0 ), ft, [&]( auto centroid, auto, auto num ) {
            for( int d = 0; d < PyPc::dim; ++d )
                ptr_res[ PyPc::dim * num + d ] = centroid[ d ];
        } );
    } );

    return res;
}

void display_vtk( const char *filename, py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid, const std::string &radial_func ) {
    VtkOutput<2> vtk_output( { "weight", "num" } );

    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    find_radial_func( radial_func, [&]( auto ft ) {
        py_grid.grid.for_each_laguerre_cell(
            [&]( auto &lc, std::size_t num_dirac_0, int ) {
                domain.bounds.for_each_intersection( lc, [&]( auto &cp, auto space_func ) {
                    cp.display( vtk_output, { ptr_weights[ num_dirac_0 ], PyPc::TF( num_dirac_0 ) } );
                } );
            }, domain.bounds.englobing_convex_polyhedron(),
            ptr_positions,
            ptr_weights,
            positions.shape( 0 ),
            false,
            need_ball_cut( ft )
        );
    } );

    vtk_output.save( filename );
}

struct PyDerResult {
    py::array_t<std::size_t> m_offsets;
    py::array_t<std::size_t> m_columns;
    py::array_t<PD_TYPE>     m_values;
    py::array_t<PD_TYPE>     v_values;
    int                      error;
};

template<class T>
void vcp( py::array_t<T> &dst, const std::vector<T> &src ) {
    dst.resize( { src.size() } );
    auto buf = dst.request();
    auto ptr = reinterpret_cast<T *>( buf.ptr );
    for( std::size_t i = 0; i < src.size(); ++i )
        ptr[ i ] = src[ i ];
}

//
PyDerResult get_der_integrals_wrt_weights( py::array_t<PD_TYPE> &positions, py::array_t<PD_TYPE> &weights, PyConvexPolyhedraAssembly &domain, PyZGrid &py_grid, const std::string &func ) {
    auto buf_positions = positions.request();
    auto buf_weights = weights.request();

    auto ptr_positions = reinterpret_cast<const PyZGrid::Pt *>( buf_positions.ptr );
    auto ptr_weights = reinterpret_cast<const PyZGrid::TF *>( buf_weights.ptr );

    std::vector<std::size_t> w_m_offsets;
    std::vector<std::size_t> w_m_columns;
    std::vector<PD_TYPE    > w_m_values;
    std::vector<PD_TYPE    > w_v_values;

    PyDerResult res;
    find_radial_func( func, [&]( auto ft ) {
        res.error = PowerDiagram::get_der_integrals_wrt_weights( w_m_offsets, w_m_columns, w_m_values, w_v_values, py_grid.grid, domain.bounds, ptr_positions, ptr_weights, std::size_t( positions.shape( 0 ) ), ft );
    } );

    vcp( res.m_offsets, w_m_offsets );
    vcp( res.m_columns, w_m_columns );
    vcp( res.m_values , w_m_values  );
    vcp( res.v_values , w_v_values  );

    return res;
}


PYBIND11_MODULE( PD_MODULE_NAME, m ) {
    m.doc() = "Power diagram tools";

    py::class_<PyZGrid>( m, "ZGrid" )
        .def( py::init<int,PD_TYPE>()                                                     , "" )
        .def( "update"                , &PyZGrid::update                                  , "" )
        .def( "display_vtk"           , &PyZGrid::display_vtk                             , "" )
    ;

    py::class_<PyConvexPolyhedraAssembly>( m, "ConvexPolyhedraAssembly" )
        .def( py::init<>()                                                                , "" )
        .def( "add_convex_polyhedron" , &PyConvexPolyhedraAssembly::add_convex_polyhedron , "" )
        .def( "add_box"               , &PyConvexPolyhedraAssembly::add_box               , "" )
        .def( "display_boundaries_vtk", &PyConvexPolyhedraAssembly::display_boundaries_vtk, "" )
        .def( "coeff_at"              , &PyConvexPolyhedraAssembly::coeff_at              , "" )
    ;    

    py::class_<PyDerResult>( m, "DerResult" )
        .def_readwrite( "m_offsets", &PyDerResult::m_offsets, "" )
        .def_readwrite( "m_columns", &PyDerResult::m_columns, "" )
        .def_readwrite( "m_values" , &PyDerResult::m_values , "" )
        .def_readwrite( "v_values" , &PyDerResult::v_values , "" )
        .def_readwrite( "error"    , &PyDerResult::error    , "" )
    ;    

    m.def( "display_vtk"                  , &display_vtk                   );
    m.def( "get_centroids"                , &get_centroids                 );
    m.def( "get_integrals"                , &get_integrals                 );
    m.def( "get_der_integrals_wrt_weights", &get_der_integrals_wrt_weights );
}

