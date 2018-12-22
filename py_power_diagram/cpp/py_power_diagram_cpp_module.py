import py_power_diagram_2d_double
import numpy as np

# import cppimport
# cppimport.force_rebuild()
# cppimport.set_quiet(False)


imported_modules = {}

# module for type and dim
def module_for_tad( type, dim ):
    if type == np.float64:
        if dim == 2 :
            return py_power_diagram_2d_double
            # name = "py_power_diagram_2d_double"
            # if not ( name in imported_modules ):
            #     imported_modules[ name ] = cppimport.imp( "py_power_diagram.internals.py_power_diagram_2d_double" )
            # return imported_modules[ name ]
        raise "TODO: bind 3D"
    raise "unmanaged type"

# module for positions and weights
def module_for_paw( positions, weights ):
    if positions.dtype != weights.dtype:
        raise "positions.dtype and weights are not of the same type"
    return module_for_tad( positions.dtype, positions.shape[ 1 ] )
