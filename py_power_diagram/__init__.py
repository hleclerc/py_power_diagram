import numpy as np
import cppimport
import os

# cppimport.force_rebuild()
# cppimport.set_quiet(False)
imported_modules = {}

def module_for( positions, weights ):
    if positions.dtype != weights.dtype:
        raise "positions.dtype and weights are not of the same type"
    if positions.dtype == np.float64:
        if positions.shape[ 1 ] == 2 :
            name = "py_power_diagram_2d_double"
            if not ( name in imported_modules ):
                imported_modules[ name ] = cppimport.imp( "py_power_diagram.py_power_diagram_2d_double" )
            return imported_modules[ name ]
        raise "TODO: bind 3D"
    raise "unmanaged positions.dtype"

class ZGrid:
    def __init__( self, max_nb_diracs_per_cell = 10 ):
        self.max_nb_diracs_per_cell = max_nb_diracs_per_cell
        self._grid = None

    def update( self, positions, weights, positions_have_changed = True, weights_have_changed = True ):
        self. _check_inst( positions, weights )
        self._grid.update( positions, weights, positions_have_changed, weights_have_changed )

    def display_vtk( self, filename ):
        if self._grid != None:
            os.makedirs( os.path.dirname( filename ), exist_ok = True )
            self._grid.display_vtk( filename )

    def _check_inst( self, positions, weights ):
        if self._grid == None:
            m = module_for( positions, weights )
            self._grid = m.ZGrid( self.max_nb_diracs_per_cell )

def grid_for( positions, weights, grid, m ):
    if grid != None:
        return grid
    res = ZGrid( 10 )
    res.update( positions, weights, True, True )
    return res

def display_vtk( positions, weights, filename, grid = None ):
    m = module_for( positions, weights )
    g = grid_for( positions, weights, grid, m )
    m.display_vtk( positions, weights, filename, g )

def get_measures( positions, weights, grid = None ):
    m = module_for( positions, weights )
    g = grid_for( positions, weights, grid, m )
    return m.get_measures( positions, weights )

