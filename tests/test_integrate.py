import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestIntegrate( unittest.TestCase ):
    def setUp( self, nb_diracs = 100 ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 1, 1 ] )

    def test_area( self, nb_diracs = 100 ):        
        for i in range( 10 ):
            # diracs
            self.positions = np.random.rand( nb_diracs, 2 )
            self.weights = np.ones( nb_diracs )

            # integrals
            areas = pd.get_integrals( "1", self.positions, self.weights, self.domain )
            self.assertAlmostEqual( np.sum( areas ), 1.0 )

if __name__ == '__main__':
    unittest.main()

# der measures
# ( has_a_void_cell, m_offsets, m_columns, m_values, v_values ) = pd.get_der_integrations_wrt_weights( positions, weights, domain )
# # display
# pd.display_vtk( "vtk/pd.vtk", positions, weights, domain )
# domain.display_boundaries_vtk( "vtk/domain.vtk" )
