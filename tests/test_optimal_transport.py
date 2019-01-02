import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestOptimalTransport( unittest.TestCase ):
    def setUp( self, nb_diracs = 100 ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 1, 1 ] )

    def test_ot( self, nb_diracs = 100 ):
        for i in range( 1 ):
            # diracs
            self.positions = np.random.rand( nb_diracs, 2 )
            self.weights = np.ones( nb_diracs )

            # optimal weights
            new_weights = pd.optimal_transport_2( "1", self.positions, self.weights, self.domain )
            print( new_weights )

            # integrals
            areas = pd.get_integrals( "1", self.positions, self.weights, self.domain )
            print( areas )
            print( np.min( areas ) )
            print( np.max( areas ) )
            # self.assertAlmostEqual( np.sum( areas ), 1.0 )

if __name__ == '__main__':
    unittest.main()
