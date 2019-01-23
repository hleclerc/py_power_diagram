import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestDerIntegrate( unittest.TestCase ):
    def setUp( self, nb_diracs = 100 ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 2, 1 ] )

    def test_unit( self ):
        self._test_der( [ 0.5, 1.5 ], [ 0, 0 ], "1" )

    def test_gaussian( self ):
        # wolfram: x=-0.5; y=(1-u)*(0.5)-0.5; N[ Integrate[ Exp[ ( 0 - x*x - y*y ) / 1 ], { u, 0, 1 } ]] -> 0.718492
        self._test_der( [ 0.4, 1.5 ], [ 1, 0 ], "exp((w-r**2)/1)" )
        self._test_der( [ 0.4, 1.5 ], [ 1, 2 ], "exp((w-r**2)/1)" )

    def _test_der( self, xs, ws, rf ):
        positions = np.array( [ [ v, 0.5 ] for v in xs ] )
        weights = np.array( [ v + 0.0 for v in ws ] )
        res = pd.get_integrals( rf, positions, weights, self.domain )
        num = pd.get_der_integrals_wrt_weights( rf, positions, weights, self.domain )

        eps = 1e-7
        for i in range( len( ws ) ):
            # self.assertAlmostEqual( num.v_values[ i ], res[ i ] )
            dweights = np.array( [ ws[ j ] + eps * ( i == j ) for j in range( len( ws ) ) ] )
            des = pd.get_integrals( rf, positions, dweights, self.domain )
            der = ( des - res ) / eps
            for j in range( len( ws ) ):
                self.assertAlmostEqual( num.m_values[ len( ws ) * i + j ], der[ j ], 6 )


if __name__ == '__main__':
    unittest.main()
