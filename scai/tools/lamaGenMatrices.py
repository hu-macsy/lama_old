#!/usr/bin/env python

import os

directory = "data/"

def generate( dim, stencil, n1, n2 = None, n3 = None ):

    filestring = "%dD%dP_%d"%( dim, stencil, n1 )
    argstring  = "%d %d %d"%(dim, stencil, n1)

    if dim > 1:
        if n2 == None:
            argstring += " %d"%n1
        else:
            filestring += "_%d"%n2
            argstring += " %d"%n2

    if dim > 2:
        if n3 == None:
            argstring += " %d"%n1
        else:
            filestring += "_%d"%n3
            argstring += " %d"%n3

    if directory:
        d = os.path.dirname( directory )
        if not os.path.exists( d ):
            os.makedirs( d )
            print 'created directory %s'%d
        filestring = d + "/" + filestring

    cmd = "./lamaGenStencilMatrix %s %s"%( filestring, argstring )

    print cmd

    if not dryrun:
        os.system( cmd )

dryrun = False

# Some useful data sets that can be used as input sets for the solvers

generate( 1, 3, 1000 )
generate( 1, 3, 1000000 )
generate( 2, 5, 500 )
generate( 2, 5, 1000 )
generate( 2, 9, 1000 )
generate( 2, 5, 2000 )
generate( 3, 7, 25 )
generate( 3, 7, 100 )
generate( 3, 19, 100 )
generate( 3, 27, 25 )
generate( 3, 27, 65 )
generate( 3, 27, 100 )
generate( 3, 27, 125 )

print ''
print ''
print 'Finished'
print '========'
print ''
print 'Input sets generated in directory %s'%directory
