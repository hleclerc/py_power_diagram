from _import_py_power_diagram import pd
import numpy as np

# domain
domain = pd.ConvexPolyhedraAssembly()
domain.add_box( [ 0, 0 ], [ 1, 1 ] )

nb_diracs = 100
positions = np.random.rand( nb_diracs, 2 )
weights   = np.ones( nb_diracs )

# display
measures = pd.get_measures( positions, weights, domain )
print( measures )
print( np.sum( measures ) )

pd.display_vtk( "vtk/pd.vtk", positions, weights, domain )
domain.display_boundaries_vtk( "vtk/domain.vtk" )
