import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np
import unittest

class TestIntegrate( unittest.TestCase ):
    def setUp( self, nb_diracs = 100 ):
        self.domain = pd.domain_types.ConvexPolyhedraAssembly()
        self.domain.add_box( [ 0, 0 ], [ 1, 1 ] )

    def test_area( self, nb_diracs = 100 ):        
        for _ in range( 10 ):
            # diracs
            positions = np.random.rand( nb_diracs, 2 )
            weights = np.ones( nb_diracs )

            # integrals
            areas = pd.get_integrals( "1", positions, weights, self.domain )
            self.assertAlmostEqual( np.sum( areas ), 1.0 )

    def test_gaussian( self ):
        # wolfram: N[ Integrate[ Integrate[ Exp[  ( 0 - x*x - y*y ) / 1 ], { x, -0.5, 0.5 } ], { y, -0.5, 0.5 } ] ]

        self._test_gaussian_for( [ 0.5, 0.5 ], w=0, eps=1.0, expected=0.851121  )
        self._test_gaussian_for( [ 0.5, 0.5 ], w=1, eps=1.0, expected=2.31359   )
        self._test_gaussian_for( [ 0.5, 0.5 ], w=0, eps=2.0, expected=0.921313  )
        self._test_gaussian_for( [ 0.5, 0.5 ], w=1, eps=2.0, expected=1.51899   )

        self._test_gaussian_for( [ 0.0, 0.0 ], w=0, eps=1.0, expected=0.557746  )
        self._test_gaussian_for( [ 0.0, 0.0 ], w=1, eps=1.0, expected=1.51611   )
        self._test_gaussian_for( [ 0.0, 0.0 ], w=0, eps=2.0, expected=0.732093  )
        self._test_gaussian_for( [ 0.0, 0.0 ], w=1, eps=2.0, expected=1.20702   )

        self._test_gaussian_for( [ 0.0, 0.0 ], w=0, eps=0.1, expected=0.0785386 )

    def _test_gaussian_for( self, positions, w, eps, expected ):
        res = pd.get_integrals( "exp((w-r**2)/{})".format( eps ), np.array( [ positions ] ), np.zeros( 1 ) + w, self.domain )
        self.assertAlmostEqual( res[ 0 ], expected, 5 )


if __name__ == '__main__':
    unittest.main()

# der measures
# ( has_a_void_cell, m_offsets, m_columns, m_values, v_values ) = pd.get_der_integrations_wrt_weights( positions, weights, domain )
# # display
# pd.display_vtk( "vtk/pd.vtk", positions, weights, domain )
# domain.display_boundaries_vtk( "vtk/domain.vtk" )
