###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief lama/solver/CMakeLists.txt
 # @author Jan Ecker
 # @date 02.04.2013
 # @since 1.0.0
###

### add subdirectories ###
add_subdirectory ( criteria )
add_subdirectory ( logger )

### for external amg ###
file ( GLOB AMG_SUBDIRPRESENT "amg" )

if ( AMG_SUBDIRPRESENT )
    add_subdirectory ( ${AMG_SUBDIRPRESENT} )
    message ( STATUS "amg subdir provides following solver libraries: ${LAMA_SOLVER_LIBRARIES}" )
    set ( LAMA_SOLVER_LIBRARIES ${LAMA_SOLVER_LIBRARIES} PARENT_SCOPE )
endif ( AMG_SUBDIRPRESENT )

### Add classes and header files ###
set ( CLASSES
        AMGSetup
        BiCG
        BiCGstab
        CG
    	CGS
    #   CGNE
        CGNR
        DefaultJacobi
        GMRES
        InverseSolver
        IterativeSolver
        MINRES
        OmegaSolver
	Richardson
        SimpleAMG
        SingleGridSetup
        SolutionProxy
        Solver
        SOR
        SpecializedJacobi
	TFQMR
        TrivialPreconditioner
    )
    
set ( HEADERS
        shared_lib
    )
    
scai_project_add ( CLASSES ${CLASSES} HEADERS ${HEADERS} SET_PARENT_SCOPE )

