import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np

import pylab

nb_diracs_par_axis = 30
target_radius = 0.15 / nb_diracs_par_axis

# domain
cp = [
     0, 0, +1, -1,
    +1, 1,  0, +1,
    -1, 1, -1, -1,
]

domain = pd.domain_types.ConvexPolyhedraAssembly()
domain.add_convex_polyhedron( cp, 1 / ( np.pi * target_radius**2 ) )
# domain.display_boundaries_vtk( "vtk/bounds.vtk" )


# initial diracs positions
def init_position_regular():
    positions = []
    for y in np.arange( 0.0, 10.0, 1.0 / nb_diracs_par_axis ):
        for x in np.arange( -10.0, 10.0, 1.0 / nb_diracs_par_axis ):
            if x**2 + y**2 <= 1 and x + y > 1.414 * target_radius and y - x > 1.414 * target_radius:
                positions.append( [ x, y ] )
    return np.array( positions )

def init_position_random():
    positions = []
    while len( positions ) < 0.5 * np.pi * nb_diracs_par_axis**2:
        x = np.random.rand() * 2 - 1
        y = np.random.rand()
        if x**2 + y**2 <= 1 and domain.contains( [ x, y ] ):
            positions.append( [ x, y ] )

    positions = np.array( positions )
    weights = np.ones( positions.shape[ 0 ] )

    loc_domain = pd.domain_types.ConvexPolyhedraAssembly()
    loc_domain.add_convex_polyhedron( cp, 1 )

    weights = pd.optimal_transport_2( "1", positions, weights, loc_domain )
    pd.display_vtk( "vtk/pd.vtk", "1", positions, weights, loc_domain )
    return pd.get_centroids( "1", positions, weights, loc_domain )



# iterations
positions = init_position_random()
weights = np.ones( positions.shape[ 0 ] ) * target_radius**2
for i in range( 100 ):
    # change positions
    positions -= 0.4 * target_radius / np.linalg.norm( positions, axis=1, keepdims=True ) * positions

    # optimal weights
    weights = pd.optimal_transport_2( "in_ball(weight**0.5)", positions, weights, domain )
    pd.display_vtk( "vtk/pd_{:03}.vtk".format( i ), "in_ball(weight**0.5)", positions, weights, domain )
  
    # update positions
    positions = pd.get_centroids( "in_ball(weight**0.5)", positions, weights, domain )

