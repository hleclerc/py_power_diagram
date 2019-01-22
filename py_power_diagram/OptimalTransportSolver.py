from petsc4py import PETSc
import numpy as np
import time

class TimeMeasurement:
    def __init__( self, table ):
        self.table = table
    def __enter__( self ):
        self.t = time.perf_counter()
    def __exit__( self, type, value, tb ):
        self.table.append( time.perf_counter() - self.t )

#
class OptimalTransportSolver:
    # 
    def __init__( self, module, domain, grid, radial_func ):
        self.radial_func = radial_func
        self.obj_max_dw  = 1e-7
        self.module      = module
        self.domain      = domain
        self.grid        = grid

        self.time_solve = []
        self.time_grid  = []
        self.time_der   = []
        self.delta_w    = []

    def solve( self, positions, weights ):
        x = PETSc.Vec().createSeq( weights.shape[ 0 ] )

        new_weights = weights + 0.0
        old_weights = weights + 0.0
        for num_iter in range( 100 ):
            # grid update
            with TimeMeasurement( self.time_grid ):
                self.grid.update( positions, new_weights, positions_have_changed = num_iter == 0, weights_have_changed = True, radial_func = self.radial_func )

            # derivatives
            with TimeMeasurement( self.time_der ):
                mvs = self.module.get_der_integrals_wrt_weights( positions, new_weights, self.domain, self.grid._inst, self.radial_func )

            # 
            if mvs.error:
                ratio = 0.1
                new_weights = ( 1 - ratio ) * old_weights + ratio * new_weights
                print( "bim (going back)" )
                continue
            old_weights = new_weights + 0.0

            #
            if self.radial_func == '1':
                mvs.m_values[ 0 ] *= 2

            with TimeMeasurement( self.time_solve ):
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

            nx = np.max( np.abs( x ) )
            new_weights -= x
            print( "max dw:", nx )

            if nx < self.obj_max_dw:
                break

        # print( 'time solve:', self.time_solve )
        # print( 'time grid :', self.time_grid  )
        # print( 'time der  :', self.time_der   )

        return new_weights


