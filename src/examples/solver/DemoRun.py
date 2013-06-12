#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Running solver example and visualization of results                    #
#                                                                         #
#    Author: Thomas Brandes, SCAI Fraunhofer                              #
#                                                                         #
#  - Measure is the runtime                                               #
#  - Also allows plotting of iterations / runtime                         #
#                                                                         #
###########################################################################

import numpy as np
import matplotlib.pyplot as plt
import subprocess

###########################################################################
#                                                                         #
#  Important configuration settings                                       #
#                                                                         #
###########################################################################

# cmd contains the parameterized command to be benchmarked
#
#   - must be a LAMA solver
#   - log level of solver must be set to advanced
#     (so output contains time for each iteration)
#
# Examples
#   cmd = "cg_solver.exe %x sp %y"
#   cmd = "cg_solver.exe 3D27P %x %y"
#
#   cg_solver.exe [inputfile] [cpu|gpu] [ell|jds|dia|csr|coo] [sync|async] [[no]texture] [no[sm]]
#

cmd     = "cg_solver.exe %x %y log_complete"

# Parameters that will be set for x and y in cmd

xlabels = ( "2D5P", "2D9P", "3D7P", "3D19P", "3D27P" )
ylabels = ( "cpu", "gpu" )
# ylabels = ( "ell", "jds", "dia" )
# ylabels = ( "32", "64", "128", "256", "512" )

# colors for different y labels

colors  = ( "r", "g", "b", "y", "m", "c", "k" )  

# enable/disable addtional plotting of time/iteration for each x value

doSinglePlot = True

# plot residual over time and not iterations

doResidualPlot = False

MaxResidual       = 10.0

# just default values, only needed for time/iteration plot

DefaultMaxTime    = 3.0
DefaultMaxIters   = 100
SampleTime        = 0.5

###########################################################################
#                                                                         #
#  Important configuration settings                                       #
#                                                                         #
###########################################################################

runningPlotFigure = None

def buildCMD( cmd, x, y ):

    # Replace %x with val of x and %y with val of y

    fullcmd = cmd.replace( "%x", x )
    fullcmd = fullcmd.replace( "%y", y )

    return fullcmd

def benchmark( cmd ):

    global MaxIters, MaxTime

    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, close_fds=True )

    niter = 0
    rtime = 0.0
    stime = 0.0

    x = list()
    y = list()

    fulltime = None
    residual = MaxResidual

    while proc.poll() is None:
        output = proc.stdout.readline()
        if "Runtime" in output:
           items = output.split( " " )
           timestr = items[ -1 ]
           fulltime = float( timestr )
           print 'Time = ', fulltime
        if "Residual" in output:
           items = output.split( " "  )
           residual =  items[-1]   # e.g. Scalar(8.8408e-06)
           residual =  residual[7:len( residual ) - 2]
           residual =  float( residual )
        if runningPlotFigure and "Duration" in output :
           items = output.split( " " )
           niter = niter + 1 
           time  = float( items[-1] )
           rtime = rtime + time
           stime = stime + time
           x.append( rtime )
           if doResidualPlot:
               y.append( residual )
           else:
               y.append( niter )
           draw = False
           if stime > SampleTime :
               if doResidualPlot:
                   plt.scatter( rtime, residual, color = plotColor )
               else:
                   plt.scatter( rtime, niter, color = plotColor )
               draw = True
               stime = 0.0
           if ( not doResidualPlot and niter > MaxIters ):
               MaxIters = 2 * MaxIters
               draw = True 
           if ( rtime > MaxTime ):
               MaxTime = 2 * MaxTime
               draw = True
           if draw:
               if doResidualPlot:
                   # plt.axis( [ 0, MaxTime, 0, MaxResidual] )
                   None
               else:
                   # plt.axis( [ 0, MaxTime, 0, MaxIters] )
                   None
               plt.draw()

    # show final point
    if runningPlotFigure:
       plt.plot( x, y, color = plotColor, linewidth = 1.0 )
       plt.draw()

       # just in case that drawing has missed the full time

       if fulltime == None:
           fulltime = rtime

    return fulltime

def measure( cmd ):

    print 'cmd = ', cmd
    time = benchmark( cmd )
    return time * 1000.0

results = []
rects   = []

title = "" 

for compVal in ylabels:
    results.append( [] )
    rects.append( [] )
    if title == "":
       title = compVal
    else:
       title = title + " / " + compVal

for x in xlabels:

    if doSinglePlot:

        runningPlotFigure = plt.figure()

        MaxTime    = DefaultMaxTime
        plt.xlabel( 'time(s)' )

        if doResidualPlot:
            MaxResidual = 10.0
            # plt.axis( [ 0, MaxTime, 0, MaxResidual] )
            plt.yscale('log')
            plt.ylabel( 'Residual' )
        else:
            MaxIters   = DefaultMaxIters
            # plt.axis( [ 0, MaxTime, 0, MaxIters] )
            plt.ylabel( '#Iterations' )

        plt.title( x + " : " + ylabels[0] + " / " + ylabels[1] )
        plt.grid( True )
        plt.ion()
        plt.show()

    for i in range( len( ylabels ) ):
   
       runCmd = buildCMD( cmd, x, ylabels[i] )

       plotColor = colors[i]
       time = measure( runCmd )
       results[i].append( time )
       print runCmd + " -> time = " + str( time )

    if doSinglePlot:
       plt.ioff()

for res in results:
   print 'Results = ', res

N = len( xlabels )
M = len( ylabels )

ind = np.arange(N)  # the x locations for the groups

width = 0.9 / len( ylabels )  # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range( len( ylabels ) ):
    rects[i] = ax.bar( ind + i * width , results[i], width, color=colors[i] )

# add some
ax.set_ylabel('Time [ms]')
ax.set_title( title )
ax.set_xticks( ind+width )
ax.set_xticklabels( xlabels )

# location for legend:  upper/lower/center left/right/center

ax.legend( rects, ylabels, loc = "upper right" )

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

for rec in rects:
    autolabel( rec )

plt.show()
