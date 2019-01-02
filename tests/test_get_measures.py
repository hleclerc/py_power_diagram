import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np

# domain
domain = pd.domain_types.ConvexPolyhedraAssembly()
domain.add_box( [ 0, 0 ], [ 1, 1 ] )

# diracs
nb_diracs = 100
positions = np.random.rand( nb_diracs, 2 )
weights   = np.ones( nb_diracs )

# integrals
areas = pd.integration( positions, weights, domain )
print( areas )
print( np.sum( areas ) )

# display
pd.display_vtk( "vtk/pd.vtk", positions, weights, domain )
domain.display_boundaries_vtk( "vtk/domain.vtk" )

# der measures
# ( has_a_void_cell, m_offsets, m_columns, m_values, v_values ) = pd.get_der_integrations_wrt_weights( positions, weights, domain )
