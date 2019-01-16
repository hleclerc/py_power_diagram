import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np

# constants
for n in [ 5 ]:
    # constants
    eps = n ** -0.5
    rfu = "exp((w-r**2)/{:.16f})".format( eps )

    # domain
    domain = pd.domain_types.ConvexPolyhedraAssembly()
    domain.add_box( [ -1, -1 ], [ 2, 2 ] )
    domain.normalize()

    domain.display_boundaries_vtk( "vtk/bounds.vtk" )

    # 
    positions = []
    for y in np.linspace( 0, 1, n ):
        for x in np.linspace( 0, 1, n ):
            positions.append( [ x, y ] )

    # iterations
    positions = np.array( positions )
    weights = np.zeros( positions.shape[ 0 ] )

    color_values = positions[ :, 1 ]
    color_values = ( color_values - np.min( color_values ) ) / ( np.max( color_values ) - np.min( color_values ) )

    #
    print( positions )

    integration = pd.get_centroids( rfu, positions, weights, domain )
    # positions = pd.get_centroids( rfu, positions, weights, domain )
    print( positions )

    # for i in range( 1 ):
    #     # change positions
    #     # positions -= 0.4 * target_radius / np.linalg.norm( positions, axis=1, keepdims=True, ord=2 ) * positions

    #     # optimal weights
    #     weights = pd.optimal_transport_2( rfu, positions, weights, domain )

    #     # display
    #     # pd.display_asy( "vtk/pd_{:03}.asy".format( i ), rfu, positions, weights, domain, values = color_values )
    #     pd.display_vtk( "vtk/pd_{:03}.vtk".format( i ), rfu, positions, weights, domain )
    
    #     # update positions
    #     positions = pd.get_centroids( rfu, positions, weights, domain )

