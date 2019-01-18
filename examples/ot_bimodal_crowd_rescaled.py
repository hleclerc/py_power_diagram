import py_power_diagram_test_context
import py_power_diagram as pd
import fast_marching as fm
import numpy as np

# constants
for na in [ 20]: # , 40, 80, 160, 200 
    directory = "vtk_{}".format( na )

    # constants
    alpha = 2*np.sqrt(1/np.pi)
    target_radius = 0.25* alpha / na
    timestep = 0.5*target_radius
    epsilon = target_radius
    
    # positions
    t = np.linspace( 0 + target_radius, alpha - target_radius, na )
    x, y = np.meshgrid( t, t )
    positions = np.hstack( ( x.reshape( ( -1, 1 ) ), y.reshape( ( -1, 1 ) ) ) )

    # domain
    domain = pd.domain_types.ConvexPolyhedraAssembly()
    domain.add_box( [ 0, 0 ], [ alpha, alpha ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.add_box( [ alpha, alpha/3 ], [ 4*alpha/3, 2*alpha/3 ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.add_box( [ 4*alpha/3, 0 ], [ 7*alpha/3, alpha ], 1 / ( np.pi * target_radius ** 2 ) )
    domain.display_boundaries_vtk( directory+"/bounds.vtk" )

    #"draw((0,0)--(3,0)--(3,1)--
    #(4,1)--(4,0)--(7,0)--(7,3)--
    #(4,3)--(4,2)--(3,2)--(3,3)--(0,3)--cycle);\n
    domain_asy = ("draw((%g,%g)--(%g,%g)--(%g,%g)--"
                  "(%g,%g)--(%g,%g)--(%g,%g)--(%g,%g)--"
                  "(%g,%g)--(%g,%g)--(%g,%g)--(%g,%g)--(%g,%g)--cycle);\n" %
                  (0, 0, alpha, 0, alpha, alpha/3,
                   4*alpha/3, alpha/3, 4*alpha/3, 0, 7*alpha/3, 0, 7*alpha/3, alpha,
                   4*alpha/3, alpha, 4*alpha/3, 2*alpha/3, alpha, 2*alpha/3, alpha, alpha, 0, alpha))
    s = 0.02
    g = fm.GradGrid( domain, [ [ 7*alpha/3-s, alpha-s ], [ 7*alpha/3-s, 0+s ] ], s )

    # iterations
    weights = np.ones( positions.shape[ 0 ] ) * target_radius ** 2
    timeout = np.zeros( positions.shape[ 0 ] )

    color_values = positions[ :, 1 ]
    color_values = ( color_values - np.min( color_values ) ) / np.ptp( color_values )
    
    h_weights = []
    h_positions = []

    T = 6.5
    nb_timesteps = int( T / timestep )
    display_timesteps = np.linspace(0, nb_timesteps-1, 100, dtype = int)
    
    for i in range( nb_timesteps ):
        print("iteration %d/%d for na=%d" % (i, nb_timesteps,na))
        
        # optimal weights
        weights = np.ones( positions.shape[ 0 ] ) * target_radius ** 2
        weights = pd.optimal_transport_2( "in_ball(weight**0.5)", positions, weights, domain )
        centroids = pd.get_centroids( "in_ball(weight**0.5)", positions, weights, domain )
        
        # display
        if i in display_timesteps:
            print(i)
            j = display_timesteps.tolist().index(i)
            pd.display_asy( directory + "/pd_{:03}.asy".format( j ), "in_ball(weight**0.5)", positions, weights, domain, values = color_values, linewidth = 0.005, dotwidth = target_radius * 0, closing = domain_asy, avoid_bounds = True )
            pd.display_vtk( directory + "/pd_{:03}.vtk".format( j ), "in_ball(weight**0.5)", positions, weights, domain )
            h_positions.append( positions )
            h_weights.append( weights )
    
        # update positions
        descent_direction = np.zeros_like(positions)
        for n in range( positions.shape[ 0 ] ):
            descent_direction[n,:] = g.grad( positions[ n, : ] )
            
        positions += timestep*(descent_direction + (centroids-positions)/epsilon)

        for n in range( positions.shape[ 0 ] ):
            if timeout[ n ] == 0 and positions[ n, 0 ] > 4:
                timeout[ n ] = i + 1

    # output with timeout information
    timeout = ( timeout - np.min( timeout ) ) / np.ptp( timeout )
    for i in range( len( h_weights ) ):
        pd.display_asy( directory + "/pd_timeout_{:03}.asy".format( i ), "in_ball(weight**0.5)", h_positions[ i ], h_weights[ i ], domain, values = timeout, linewidth = 0.005, dotwidth = target_radius * 0, closing = domain_asy, avoid_bounds = True )
