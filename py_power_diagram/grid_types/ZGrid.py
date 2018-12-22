from py_power_diagram.cpp import py_power_diagram_cpp_module
import numpy as np
import os

# grid based on morton ordering
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
            m = py_power_diagram_cpp_module.module_for_paw( positions, weights )
            self._inst = m.ZGrid( self.max_nb_diracs_per_cell )
