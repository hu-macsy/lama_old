#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Benchmarking and Visualization                                         #
#                                                                         #
###########################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import subprocess

###########################################################################
#                                                                         #
#  Global variables use                                                   #
#                                                                         #
###########################################################################

# colors for different y labels

colors  = ( "r", "g", "b", "y", "m", "c", "k" )  

# enable/disable addtional plotting of time/iteration for each x value

doSinglePlot = True

# plot residual over time and not iterations

doResidualPlot = True

# interval time: after each this time a new point is written in the single plots

SampleTime        = 2.0

# String that must appear in an out, e.g. for MPI: MPI(0:1), to select node

selectOutput      = None

###########################################################################
#                                                                         #
#  Interface                                                              #
#                                                                         #
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

# xlist contains values to be set for %s in cmd
# ylist contains values to be set for %s in cmd

# singleRunChoice = 0    no individual visualization for each x values

def evalSingleRunChoice( singleRunChoice ):

    global doSinglePlot, doResidualPlot

    if singleRunChoice == 0:

        doSinglePlot = False

    elif singleRunChoice == 1:

        doSinglePlot = True
        doResidualPlot = False
       
    elif singleRunChoice == 2:

        doSinglePlot = True
        doResidualPlot = True

    else:

        print "Illegal argument for singleRunChoice (0 = none, 1 = Iters, 2 = Residual)"
        raise ArgumentError( singleRunChoice )

def buildCMD( cmd, x, y ):

    # Replace all %x with val of x and all %y with val of y

    cmdX = cmd.replace( "%x", str( x ) )

    if cmdX == cmd:
        print "Attention: cmd = %s does not contain %s"%( cmd, "%x" )

    cmdXY = cmdX.replace( "%y", str( y ) )

    if cmdXY == cmdX:
        print "Attention: cmd = %s does not contain %s"%( cmd, "%y" )

    return cmdXY

niter = 0
rtime = 0.0
stime = 0.0

def evalOutput( output ):

    global niter, fulltime, rtime, stime, x, y

    if outputSelection != None:

       if not outputSelection in output:
          return

    if "Runtime" in output:

       items = output.split( " " )
       timestr = items[ -1 ]
       fulltime = float( timestr )
       print 'Runtime = ', fulltime

    if not doSinglePlot:
       return 0.0

    if "Residual" in output:
       items = output.split( " "  )
       residual =  items[-1]   # e.g. 2.221\n
       residual =  float( residual )
       if doResidualPlot:
          x.append( rtime )
          y.append( residual )

    if "Duration" in output :
       items = output.split( " " )
       niter = niter + 1 
       time  = float( items[-1] )
       rtime = rtime + time
       stime = stime + time
       if not doResidualPlot:
          x.append( rtime )
          y.append( niter )

    if stime > SampleTime :

       # scatter last value to see progress

       if len( x ) > 0:
           plt.scatter( x[-1], y[-1], color = plotColor )
       plt.pause(0.0001)
       stime = 0.0

def benchmark( cmd ):

    global fulltime
    global x, y, niter, rtime, stime

    fulltime = None
    niter = 0
    rtime = 0.0
    stime = 0.0
    
    x = list()
    y = list()

    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, close_fds=True )

    while proc.poll() is None:

        # evaluate output lines of busy process

        output = proc.stdout.readline()
        evalOutput( output )

    # There might be some missing lines at the end of the process

    for output in proc.stdout:
        evalOutput( output )

    # show final point
    if doSinglePlot:
       plt.plot( x, y, color = plotColor, linewidth = 1.0 )
       plt.pause(0.0001)

       # just in case that drawing has missed the full time

       if fulltime == None:
           fulltime = rtime

    return fulltime

def measure( cmd ):

    print 'cmd = ', cmd
    time = benchmark( cmd )
    return time * 1000.0

def autolabel( ax, rects ):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text( rect.get_x() + rect.get_width()/2., 1.05*height, '%d'%int(height),
                 ha='center', va='bottom')

def setupSinglePlot( title ):
    
    # set up a new plot figure for  time / residual or  time / #iterations

    plt.xlabel( 'time(s)' )
    
    if doResidualPlot:
        plt.yscale( 'log' )
        plt.ylabel( 'Residual' )
    else:
        plt.ylabel( '#Iterations' )
        plt.axis( [ 0, 1, 0, 10 ] )
        ax = plt.gca()
        ax.set_autoscale_on( True )

    plt.title( title )
    plt.grid( True )
    plt.ion()
    plt.show()

plotColor = 'g'

def MultiBench( cmd, xlist, ylist, singleRunChoice, title = None, selection = None, save = None ):

    global plotColor
    global outputSelection

    # make some global settings for single runs

    evalSingleRunChoice( singleRunChoice )

    outputSelection = selection

    fig = plt.figure( figsize = ( 24.0, 16.0 ) )

    fig.canvas.set_window_title( 'LAMA : ' + str( title ) )

    plt.suptitle( str( title ) )

    plotid = 1
    plotRows = 1
    plotColumns = 1
    plotN       = 1

    if doSinglePlot:
       plotN = 1 + len( xlist )
       plotColumns = plotN

    if plotColumns > 2:
       plotRows = 2
       plotColumns = ( plotColumns + 1 ) / 2

    print "Plot grid = %d x %d, plots = %d"%( plotRows, plotColumns, plotN )

    lamaImage = mpimg.imread( "lama.png" )

    results = []
    rects   = []

    title = "" 

    for compVal in ylist:

        results.append( [] )
        rects.append( [] )

        if title == "":
           title = str( compVal )
        else:
           title = title + " / " + str( compVal )

    for x in xlist:
    
        if doSinglePlot:
   
            fig.add_subplot( plotRows, plotColumns, plotid )
            plotid = plotid + 1 
            setupSinglePlot( x + " : " + title )
    
        for i in range( len( ylist ) ):
       
            runCmd = buildCMD( cmd, x, ylist[i] )
    
            plotColor = colors[i]
            time = measure( runCmd )
            results[i].append( time )
            print runCmd + " -> time = " + str( time )
    
        if doSinglePlot:

            plt.ioff()
    
    for res in results:
        print 'Results = ', res

    N = len( xlist )
    M = len( ylist )

    ind = np.arange(N)  # the x locations for the groups

    width = 0.9 / len( ylist )  # the width of the bars
    
    ax = fig.add_subplot( plotRows, plotColumns, plotid )
    plotid = plotid + 1

    for i in range( len( ylist ) ):
        rects[i] = ax.bar( ind + i * width , results[i], width, color=colors[i] )

    # add some
    ax.set_ylabel('Time [ms]')
    ax.set_title( title )
    ax.set_xticks( ind + width )
    ax.set_xticklabels( xlist )

    # location for legend:  upper/lower/center left/right/center

    ax.legend( rects, ylist, loc = "upper right" )

    for rec in rects:
        autolabel( ax, rec )
    
    if ( plotRows * plotColumns ) != plotN:
       fig.add_subplot( plotRows, plotColumns, plotid )
       plotid = plotid + 1 
       plt.imshow( lamaImage )

    if save != None:
       plt.savefig( save )
       print "Saved figure as %s"%save

    plt.savefig( "lastPic.png", bbox_inches = 0 )

    plt.show()
