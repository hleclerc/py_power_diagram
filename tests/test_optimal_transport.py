import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestOptimalTransport( unittest.TestCase ):
    def setUp( self ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 1, 1 ] )

    def test_ot( self, nb_diracs = 100000 ):
        for _ in range( 1 ):
            # diracs
            self.positions = np.random.rand( nb_diracs, 2 )
            self.weights = np.ones( nb_diracs )

            # optimal weights
            new_weights = pd.optimal_transport_2( "1", self.positions, self.weights, self.domain )

            # integrals
            areas = pd.get_integrals( "1", self.positions, new_weights, self.domain )
            self.assertAlmostEqual( np.min( areas ), 1.0 / nb_diracs, places = 6 )
            self.assertAlmostEqual( np.max( areas ), 1.0 / nb_diracs, places = 6 )

            # pd.display_vtk( "vtk/pd.vtk", self.positions, self.weights, self.domain )

if __name__ == '__main__':
    np.random.seed( 1 )
    unittest.main()
