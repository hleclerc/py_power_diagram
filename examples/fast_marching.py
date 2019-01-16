import numpy as np
import skfmm
import math

import matplotlib.pyplot

class GradGrid:
    def __init__( self, domain, targets, space_step ):
        self.mi = domain.min_position() - space_step
        self.ma = domain.max_position() + space_step
        self.ss = space_step

        nd = self.mi.shape[ 0 ]

        # make a grid
        n = [ math.ceil( ( self.ma[ d ] - self.mi[ d ] ) / space_step )  for d in range( nd ) ]
        s = [ np.linspace( self.mi[ d ], self.ma[ d ], n[ d ] ) for d in range( nd ) ]
        g = np.meshgrid( *s, indexing='ij' )
        c = np.hstack( [ g[ d ].reshape( ( -1, 1 ) ) for d in range( nd ) ] )

        # forbidden zones
        mask = np.array( [ not domain.contains( c[ i, : ] ) for i in range( c.shape[ 0 ] ) ] ).reshape( n )
        
        # forbidden zones
        phi = - 1.0 * np.ones_like( mask ) # 
        for t in targets:
            c = np.floor( ( t - self.mi ) / space_step ).astype( int )
            if nd == 3:
                raise "TODO"
            phi[ c[ 0 ], c[ 1 ] ] = 0

        self.dist = skfmm.distance( np.ma.MaskedArray( phi, mask ), dx = [ space_step for d in range( nd ) ] )

        # precomputation of an integer spiral (used in grad func)
        self._offsets = []
        if nd == 2:
            for y in range( -5, 6 ):
                for x in range( -5, 6 ):
                    self._offsets.append( ( x, y ) )
            self._offsets.sort( key = lambda p: p[ 0 ]**2 + p[ 1 ]**2 )
        else:
            raise "TODO"

    def grad( self, position ):
        ib = np.floor( ( position - self.mi ) / self.ss ).astype( int )
        nd = self.mi.shape[ 0 ]
        if nd == 3:
            raise "TODO"
        for dx, dy in self._offsets:
            ix = ib[ 0 ] + dx
            iy = ib[ 1 ] + dy
            if ix >= 0 and iy >= 0 and ix < np.shape( self.dist.mask )[ 0 ] - 1 and iy < np.shape( self.dist.mask )[ 1 ] - 1 and self.dist.mask[ ix, iy ] == False and self.dist.mask[ ix + 1, iy ] == False and self.dist.mask[ ix, iy + 1 ] == False:
                gx = self.dist[ ix + 1, iy ] - self.dist[ ix, iy ]
                gy = self.dist[ ix, iy + 1 ] - self.dist[ ix, iy ]
                re = np.array( [ gx, gy ] )
                return re / np.linalg.norm( re )
        return np.array( [ 0, 1 ] )
