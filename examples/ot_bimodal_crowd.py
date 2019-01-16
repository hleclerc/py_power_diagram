import py_power_diagram_test_context
import py_power_diagram as pd
import fast_marching as fm
import numpy as np

# constants
for n in [ 30 ]:
    # constants
    target_radius = 3 * 0.45 / n

    # positions
    t = np.linspace( 0 + target_radius, 3 - target_radius, n )
    x, y = np.meshgrid( t, t )
    positions = np.hstack( ( x.reshape( ( -1, 1 ) ), y.reshape( ( -1, 1 ) ) ) )

    # domain
    domain = pd.domain_types.ConvexPolyhedraAssembly()
    domain.add_box( [ 0, 0 ], [ 3, 3 ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.add_box( [ 3, 1 ], [ 4, 2 ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.add_box( [ 4, 0 ], [ 7, 3 ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.display_boundaries_vtk( "vtk/bounds.vtk" )

    s = 0.02
    g = fm.GradGrid( domain, [ [ 7-s, 3-s ], [ 7-s, 0+s ] ], s )

    # iterations
    weights = np.ones( positions.shape[ 0 ] ) * target_radius ** 2

    color_values = positions[ :, 1 ]
    color_values = ( color_values - np.min( color_values ) ) / ( np.max( color_values ) - np.min( color_values ) )

    for i in range( int( 20 / target_radius ) ):
        # change positions
        for n in range( positions.shape[ 0 ] ):
            positions[ n, : ] += 0.5 * target_radius * g.grad( positions[ n, : ] )

        # optimal weights
        weights = pd.optimal_transport_2( "in_ball(weight**0.5)", positions, weights, domain )

        # display
        pd.display_asy( "vtk/pd_{:03}.asy".format( i ), "in_ball(weight**0.5)", positions, weights, domain, values = color_values, linewidth = 0.01, dotwidth = target_radius * 0 )
        pd.display_vtk( "vtk/pd_{:03}.vtk".format( i ), "in_ball(weight**0.5)", positions, weights, domain )
    
        # update positions
        positions = pd.get_centroids( "in_ball(weight**0.5)", positions, weights, domain )
