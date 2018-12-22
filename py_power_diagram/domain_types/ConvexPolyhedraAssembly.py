from py_power_diagram.cpp import py_power_diagram_cpp_module
import numpy as np
import os

#
class ConvexPolyhedraAssembly:
    def __init__( self, dim = 2, type = np.float64 ):
        self._inst = py_power_diagram_cpp_module.module_for_tad( type, dim ).ConvexPolyhedraAssembly()
        self._type = type

    def add_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self._inst.add_box( np.array( min_pos, dtype=self._type ), np.array( max_pos, dtype=self._type ), self._type( coeff ), np.uint64( cut_id ) )

    def sub_box( self, min_pos, max_pos, coeff = 1.0, cut_id = -1 ):
        self.add_box( min_pos, max_pos, - coeff, cut_id )

    def display_boundaries_vtk( self, filename ):
        os.makedirs( os.path.dirname( filename ), exist_ok = True )
        self._inst.display_boundaries_vtk( filename )