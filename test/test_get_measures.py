from _import_py_power_diagram import pd
import numpy as np

nb_diracs = 10
positions = np.random.rand( nb_diracs, 2 )
weights   = np.ones( nb_diracs )

pd.display_vtk( positions, weights, "vtk/pd.vtk" )

