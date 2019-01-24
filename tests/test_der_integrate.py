import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestDerIntegrate( unittest.TestCase ):
    def setUp( self, nb_diracs = 100 ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 2, 1 ] )

    def test_unit( self ):
        self._test_der( [ [ 0.5, 0.5 ], [ 1.5, 0.5 ] ], [ 0, 0 ], "1" )

    def test_gaussian( self ):
        # wolfram: x=-0.5; y=(1-u)*(0.5)-0.5; N[ Integrate[ Exp[ ( 0 - x*x - y*y ) / 1 ], { u, 0, 1 } ]] -> 0.718492
        self._test_der( [ [ 0.4, 0.5 ], [ 1.5, 0.5 ] ], [ 1, 0 ], "exp((w-r**2)/1)" )
        self._test_der( [ [ 0.4, 0.5 ], [ 1.5, 0.5 ] ], [ 1, 2 ], "exp((w-r**2)/1)" )


        N = 10
        for i in range( 50 ):
            positions = [ [ np.random.rand(), np.random.rand() ] for i in range( N ) ]
            weights = [ np.random.rand() * 0.1 for i in range( N ) ]
            self._test_der( positions, weights, "exp((w-r**2)/2)" )


    def _test_der( self, positions, weights, rf ):
        weights = np.array( weights, dtype=np.double )
        positions = np.array( positions, dtype=np.double )
        num = pd.get_der_integrals_wrt_weights( rf, positions, weights, self.domain )
        if num.error:
            return

        mat = np.zeros( ( positions.shape[ 0 ], positions.shape[ 0 ] ) )
        for i in range( positions.shape[ 0 ] ):
            for o in range( num.m_offsets[ i + 0 ], num.m_offsets[ i + 1 ] ):
                mat[ i, num.m_columns[ o ] ] = num.m_values[ o ]

        eps = 1e-6
        delta = np.max( np.abs( mat ) ) * 100 * eps
        res = pd.get_integrals( rf, positions, weights, self.domain )
        for i in range( positions.shape[ 0 ] ):
            dweights = np.array( [ weights[ j ] + eps * ( i == j ) for j in range( positions.shape[ 0 ] ) ] )
            des = pd.get_integrals( rf, positions, dweights, self.domain )
            der = ( des - res ) / eps
            for j in range( positions.shape[ 0 ] ):
                self.assertAlmostEqual( mat[ i, j ], der[ j ], delta=delta )


if __name__ == '__main__':
    unittest.main()
