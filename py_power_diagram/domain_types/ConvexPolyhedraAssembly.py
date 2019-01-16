from py_power_diagram.cpp import py_power_diagram_cpp_module
import numpy as np
import os

#
class ConvexPolyhedraAssembly:
    def __init__( self, dim = 2, type = np.float64 ):
        self._inst = py_power_diagram_cpp_module.module_for_tad( type, dim ).ConvexPolyhedraAssembly()
        self._type = type
        self._dim  = dim

    def add_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self._inst.add_box( np.array( min_pos, dtype=self._type ), np.array( max_pos, dtype=self._type ), self._type( coeff ), np.uint64( cut_id ) )

    def sub_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self.add_box( min_pos, max_pos, - coeff, cut_id )

    # remember to call normalize when integration( coeff ) != 1
    def add_convex_polyhedron( self, positions_and_normals, coeff = 1.0, cut_id = -1 ):
        pan = np.array( positions_and_normals, dtype=self._type )
        if len( pan.shape ) == 1:
            pan = pan.reshape( -1, 2 * self._dim )
        self._inst.add_convex_polyhedron( pan, self._type( coeff ), np.uint64( cut_id ) )

    def normalize( self ):
        self._inst.normalize()

    def display_boundaries_vtk( self, filename ):
        os.makedirs( os.path.dirname( filename ), exist_ok = True )
        self._inst.display_boundaries_vtk( filename )

    # coefficient at `point`. If point is not contained, return 0.
    def coeff_at( self, point ):
        return self._inst.coeff_at( np.array( point, dtype=self._type ) )

    # 
    def min_position( self ):
        return self._inst.min_position()

    # 
    def max_position( self ):
        return self._inst.max_position()

    def measure( self ):
        return self._inst.measure()

    # True if points is contained
    def contains( self, point ):
        return self.coeff_at( point ) != 0
