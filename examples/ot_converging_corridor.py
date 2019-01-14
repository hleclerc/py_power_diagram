import py_power_diagram_test_context
import py_power_diagram as pd
import numpy as np

import pylab

nb_diracs_par_axis = 10
target_radius = 0.1

# domain
domain = pd.domain_types.ConvexPolyhedraAssembly()
domain.add_convex_polyhedron( [
     0, 0, +1, -1,
    +9, 9,  0, +1,
    -9, 9, -1, -1,
], 1 / ( 3.14159 * target_radius**2 ) )


# diracs
positions = []
weights = []
for y in np.arange( 0.0, 10.0, 1.0 / nb_diracs_par_axis ):
    for x in np.arange( -10.0, 10.0, 1.0 / nb_diracs_par_axis ):
        if x**2 + y**2 <= 1 and domain.contains( [ x, y ] ):
            positions.append( [ x, y ] )
            weights.append( 0.5 )
positions = np.array( positions )
weights = np.array( weights )

# optimal weights
# new_weights = pd.optimal_transport_2( "in_ball(weight**0.5)", positions, weights, domain )
# new_weights = pd.optimal_transport_2( "1", positions, weights, domain )

# domain.display_boundaries_vtk( "vtk/bounds.vtk" )
# pd.display_vtk( "vtk/pd.vtk", "in_ball(weight**0.5)", positions, new_weights, domain )
pd.display_vtk( "vtk/pd.vtk", "in_ball(weight**0.5)", positions, weights, domain )
