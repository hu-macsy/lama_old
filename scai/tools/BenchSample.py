#!/usr/bin/env python

import MultiBench

# Important: log_complete must be switched on for Logger fo Solver

cmd     = "./lamaSolver data/3D27P_100.frm --SCAI_FORMAT=%x --SCAI_NUM_THREADS=%y --SCAI_SOLVER_LOG=completeInformation"

xlabels = ( "CSR", "ELL", "JDS", "DIA" )
ylabels = ( "1", "2" )

# Just final results

# kind = 0

# One figure for each x label, plot iterations

# kind = 1

# One figure for each x label, plot residual

kind = 2

title = "CG solver of 2D5P stencil 500 x 500 - different sparse matrix formats"

MultiBench.MultiBench( cmd, xlabels, ylabels, kind, title = title, save = "BenchFormat2D.png" )

