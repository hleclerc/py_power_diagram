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
        positions = np.array( [[ 0.5, 0.5 ]] )

        # wolfram: N[ Integrate[ Integrate[ Exp[  ( 0 - x*x - y*y ) / 1 ], { x, -0.5, 0.5 } ], { y, -0.5, 0.5 } ] ]
        res = pd.get_integrals( "exp((w-r**2)/1)", positions, np.zeros( 1 ), self.domain )
        self.assertAlmostEqual( res[ 0 ], 0.851121, 5 )

        # wolfram: N[ Integrate[ Integrate[ Exp[  ( 1 - x*x - y*y ) / 1 ], { x, -0.5, 0.5 } ], { y, -0.5, 0.5 } ] ]
        res = pd.get_integrals( "exp((w-r**2)/1)", positions, np.ones( 1 ), self.domain )
        self.assertAlmostEqual( res[ 0 ], 2.31359, 5 )

        # # wolfram: N[ Integrate[ Integrate[ Exp[  ( 1 - x*x - y*y ) / 1 ], { x, -0.5, 0.5 } ], { y, -0.5, 0.5 } ] ]
        # res = pd.get_integrals( "exp((w-r**2)/1)", positions, np.ones( 1 ), self.domain )
        # self.assertAlmostEqual( res[ 0 ], 2.31359, 5 )
        
        # same thing with a different position
        positions = np.array( [[ 0.0, 0.0 ]] )

        # wolfram: N[ Integrate[ Integrate[ Exp[  ( 0 - x*x - y*y ) / 1 ], { x, 0, 1 } ], { y, 0, 1 } ] ]
        res = pd.get_integrals( "exp((w-r**2)/1)", positions, np.zeros( 1 ), self.domain )

        print( "integration:", res )

if __name__ == '__main__':
    unittest.main()

# der measures
# ( has_a_void_cell, m_offsets, m_columns, m_values, v_values ) = pd.get_der_integrations_wrt_weights( positions, weights, domain )
# # display
# pd.display_vtk( "vtk/pd.vtk", positions, weights, domain )
# domain.display_boundaries_vtk( "vtk/domain.vtk" )
