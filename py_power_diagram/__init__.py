from .OptimalTransportSolver import OptimalTransportSolver
from . import domain_types
from . import grid_types
from . import cpp
import numpy as np
import os

def _grid_for( positions, weights, grid, m ):
    if grid != None:
        if grid._inst == None:
            grid.update( positions, weights )
        return grid
    res = grid_types.ZGrid( 13 )
    res.update( positions, weights )
    return res

def _domain_for( positions, weights, domain, m ):
    if domain != None:
        return domain._inst
    res = domain_types.ConvexPolyhedraAssembly()
    res.add_box( [ 0, 0 ], [ 1.5, 1 ], 1.0, -1 )
    return res._inst

# make a vtk file for a representation of the power diagram
def display_vtk( filename, func, positions, weights, domain = None, grid = None ):
    os.makedirs( os.path.dirname( filename ), exist_ok = True )
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    m.display_vtk( filename, positions, weights, d, g._inst, func.lower() )

# make a metapost file for a representation of the power diagram
def display_asy( filename, func, positions, weights, domain = None, grid = None, preamble = "", output_format = "pdf", linewidth = 0.02, dotwidth = 0.0, values = np.array([]), colormap = "inferno" ):
    os.makedirs( os.path.dirname( filename ), exist_ok = True )
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )

    p = "settings.outformat = \"{}\";\nunitsize(1cm);\n".format( output_format )
    if linewidth > 0:
        p += "defaultpen({}cm);\n".format( linewidth )
    elif dotwidth > 0:
        p += "defaultpen({}cm);\n".format( dotwidth / 6 )
    p += preamble

    m.display_asy( filename, positions, weights, d, g._inst, func.lower(), p, values, colormap, linewidth, dotwidth )

# return integral( cell_i, func... ) for each cell
# possible func types: "1", ...
def get_integrals( func, positions, weights, domain = None, grid = None ):
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    return m.get_integrals( positions, weights, d, g._inst, func.lower() )

# return centroid( cell_i, func... ) for each cell
# possible func types: "1", ...
def get_centroids( func, positions, weights, domain = None, grid = None ):
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    return m.get_centroids( positions, weights, d, g._inst, func.lower() )

# derivatives of integrals of $func wrt weights
# possible func types: "1", ...
def get_der_integrals_wrt_weights( func, positions, weights, domain = None, grid = None ):
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    return m.get_der_integrals_wrt_weights( positions, weights, d, g._inst, func.lower() )

# return optimal weights
# possible func types: "1", ...
def optimal_transport_2( func, positions, weights, domain = None, grid = None ):
    m = cpp.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )

    s = OptimalTransportSolver( m, d, g, func.lower() )
    return s.solve( positions, weights )


