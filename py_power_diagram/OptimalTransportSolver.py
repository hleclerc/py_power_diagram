from petsc4py import PETSc
import numpy as np
import time

#
class OptimalTransportSolver:
    def __init__( self, module, domain, grid, func ):
        self.obj_max_dw = 1e-7
        self.module     = module
        self.domain     = domain
        self.grid       = grid
        self.func       = func

        self.time_solve = []
        self.time_grid  = []
        self.time_der   = []
        self.delta_w    = []

    def solve( self, positions, weights ):
        x = PETSc.Vec().createSeq( weights.shape[ 0 ] )

        new_weights = weights + 0.0
        old_weights = weights + 0.0
        for num_iter in range( 100 ):
            self.delta_w.append( np.max( new_weights ) - np.min( new_weights ) )

            # grid update
            t0 = time.perf_counter()
            self.grid.update( positions, new_weights, num_iter == 0, True )

            t1 = time.perf_counter()
            self.time_grid.append( t1 - t0 )

            # derivatives
            mvs = self.module.get_der_integrals_wrt_weights( positions, new_weights, self.domain, self.grid, self.func )

            t2 = time.perf_counter()
            self.time_der.append( t2 - t1 )

            # 
            if mvs.error:
                ratio = 0.1
                new_weights = ( 1 - ratio ) * old_weights + ratio * new_weights
                print( "bim" )
                continue
            old_weights = new_weights + 0.0

            #
            if self.func == '1':
                mvs.m_values[ 0 ] *= 2

            A = PETSc.Mat().createAIJ( [ weights.shape[ 0 ], weights.shape[ 0 ] ], csr = ( mvs.m_offsets.astype( PETSc.IntType ), mvs.m_columns.astype( PETSc.IntType ), mvs.m_values ) )
            b = PETSc.Vec().createWithArray( mvs.v_values ) 
            A.assemblyBegin() # Make matrices useable.
            A.assemblyEnd()
            
            # Initialize ksp solver.
            ksp = PETSc.KSP().create()
            ksp.setType( 'cg' )
            # ksp.getPC().setType( 'icc' )
            ksp.getPC().setType( 'gamg' )
            # print( 'Solving with:', ksp.getType() )

            ksp.setOperators( A )
            ksp.setFromOptions()
            
            # Solve
            ksp.solve( b, x )

            t3 = time.perf_counter()
            self.time_solve.append( t3 - t2 )

            nx = np.max( np.abs( x ) )
            new_weights -= x
            print( nx )

            if nx < self.obj_max_dw:
                break

        print( 'solve:', self.time_solve )
        print( 'grid :', self.time_grid  )
        print( 'der  :', self.time_der   )
        print( 'dw   :', self.delta_w    )

        print( 'solve:', np.sum( self.time_solve ) )
        print( 'grid :', np.sum( self.time_grid  ) )
        print( 'der  :', np.sum( self.time_der   ) )
        return new_weights


