#!/usr/bin/env python

import MultiBench

import os

matrix1 = "matrix_generator.exe 2D9P 2 9 500 500"
matrix2 = "matrix_generator.exe 3D27P 3 27 65 65 65"

os.system( matrix1 )
os.system( matrix2 )

print 'Input sets generated'

# Important: log_complete must be set

cmd     = "cg_solver.exe %x cpu csr %y cpu log_complete gpu"

xlabels = ( "2D9P", "3D27P" )
ylabels = ( "csr", "dia", "ell", "jds" )

# Just final results
# MultiBench.MultiBench( cmd, xlabels, ylabels, 0 )

# One figure for each x label, plot iterations
# MultiBench.MultiBench( cmd, xlabels, ylabels, 1 )

# One figure for each x label, plot residual
MultiBench.MultiBench( cmd, xlabels, ylabels, 2 )
