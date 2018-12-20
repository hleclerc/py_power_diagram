import numpy as np
import cppimport
import os

# cppimport.force_rebuild()
# cppimport.set_quiet(False)
imported_modules = {}

def _module_for_tad( type, dim ):
    if type == np.float64:
        if dim == 2 :
            name = "py_power_diagram_2d_double"
            if not ( name in imported_modules ):
                imported_modules[ name ] = cppimport.imp( "py_power_diagram.py_power_diagram_2d_double" )
            return imported_modules[ name ]
        raise "TODO: bind 3D"
    raise "unmanaged type"

def _module_for_paw( positions, weights ):
    if positions.dtype != weights.dtype:
        raise "positions.dtype and weights are not of the same type"
    return _module_for_tad( positions.dtype, positions.shape[ 1 ] )

class ConvexPolyhedraAssembly:
    def __init__( self, dim = 2, type = np.float64 ):
        self._inst = _module_for_tad( type, dim ).ConvexPolyhedraAssembly()
        self._type = type

    def add_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self._inst.add_box( np.array( min_pos, dtype=self._type ), np.array( max_pos, dtype=self._type ), self._type( coeff ), np.uint64( cut_id ) )

    def sub_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self.add_box( min_pos, max_pos, - coeff, cut_id )

    def display_boundaries_vtk( self, filename ):
        os.makedirs( os.path.dirname( filename ), exist_ok = True )
        self._inst.display_boundaries_vtk( filename )

class ZGrid:
    def __init__( self, max_nb_diracs_per_cell = 10 ):
        self.max_nb_diracs_per_cell = max_nb_diracs_per_cell
        self._inst = None

    def update( self, positions, weights, positions_have_changed = True, weights_have_changed = True ):
        self. _check_inst( positions, weights )
        self._inst.update( positions, weights, positions_have_changed, weights_have_changed )

    def display_vtk( self, filename ):
        if self._inst != None:
            os.makedirs( os.path.dirname( filename ), exist_ok = True )
            self._inst.display_vtk( filename )

    def _check_inst( self, positions, weights ):
        if self._inst == None:
            m = _module_for_paw( positions, weights )
            self._inst = m.ZGrid( self.max_nb_diracs_per_cell )

def _grid_for( positions, weights, grid, m ):
    if grid != None:
        return grid._inst
    res = ZGrid( 10 )
    res.update( positions, weights, True, True )
    return res._inst

def _domain_for( positions, weights, domain, m ):
    if domain != None:
        return domain._inst
    res = ConvexPolyhedraAssembly()
    res.add_box( [ 0, 0 ], [ 1.5, 1 ], 1.0, -1 )
    return res._inst

def display_vtk( filename, positions, weights, domain = None, grid = None ):
    os.makedirs( os.path.dirname( filename ), exist_ok = True )
    m = _module_for_paw( positions, weights )
    g = _grid_for( positions, weights, grid, m )
    d = _domain_for( positions, weights, domain, m )
    m.display_vtk( filename, positions, weights, d, g )

def get_measures( positions, weights, domain = None, grid = None ):
    m = _module_for_paw( positions, weights )
    g = _grid_for( positions, weights, grid, m )
    d = _domain_for( positions, weights, domain, m )
    return m.get_measures( positions, weights, d, g )

