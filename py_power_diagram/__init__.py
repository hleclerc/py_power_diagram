from . import domain_types
from . import grid_types
from . import internals
import numpy as np
import os

def _grid_for( positions, weights, grid, m ):
    if grid != None:
        if grid._inst == None:
            grid.update( positions, weights )
        return grid._inst
    res = grid_types.ZGrid( 10 )
    res.update( positions, weights )
    return res._inst

def _domain_for( positions, weights, domain, m ):
    if domain != None:
        return domain._inst
    res = domain_types.ConvexPolyhedraAssembly()
    res.add_box( [ 0, 0 ], [ 1.5, 1 ], 1.0, -1 )
    return res._inst

# make a vtk file for a representation of the power diagram
def display_vtk( filename, positions, weights, domain = None, grid = None ):
    os.makedirs( os.path.dirname( filename ), exist_ok = True )
    m = internals.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    m.display_vtk( filename, positions, weights, d, g )

# return integral( cell_i d... ) for each cell
# possible measure types: "leb", ...
def integration( positions, weights, domain = None, grid = None, func = "1" ):
    m = internals.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    return m.integration( positions, weights, d, g, func.lower() )

# wrt weights
def der_integration_wrt_weights( positions, weights, domain = None, grid = None, func = "1" ):
    m = internals.py_power_diagram_cpp_module.module_for_paw( positions, weights )
    d = _domain_for( positions, weights, domain, m )
    g = _grid_for( positions, weights, grid, m )
    return m.get_der_integrations( positions, weights, d, g, func.lower() )

