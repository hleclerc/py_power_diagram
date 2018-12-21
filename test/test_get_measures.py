from _import_py_power_diagram import pd
import numpy as np

# domain
domain = pd.domain_types.ConvexPolyhedraAssembly()
domain.add_box( [ 0, 0 ], [ 1, 1 ] )

np.random.seed( 5 )

# diracs
nb_diracs = 100
positions = np.random.rand( nb_diracs, 2 )
weights   = np.ones( nb_diracs )

grid = pd.grid_types.ZGrid( 1000 )
grid.update( positions, weights )
grid.display_vtk( "vtk/grid.vtk" )

# integrals
areas = pd.integration( positions, weights, domain, grid )
print( areas )
print( np.sum( areas ) )

# display
pd.display_vtk( "vtk/pd.vtk", positions, weights, domain, grid )
domain.display_boundaries_vtk( "vtk/domain.vtk" )

# der measures
# ( has_a_void_cell, m_offsets, m_columns, m_values, v_values ) = pd.get_der_integrations_wrt_weights( positions, weights, domain )
