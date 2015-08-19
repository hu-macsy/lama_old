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
 # @brief Everything needed for using Graph partitioned Distribution ( only with Metis/ParMetis yet )
 # @author Lauretta Schubert
 # @date 19.08.2015
 # @since 2.0.0
###

find_package ( Metis ${LAMA_FIND_PACKAGE_FLAGS} )
if    ( METIS_FOUND )
	find_package ( ParMetis ${LAMA_FIND_PACKAGE_FLAGS} )
endif ( METIS_FOUND )

<<<<<<< HEAD:CMake/Modules/Package/GraphPartitioning.cmake
## ALLOW to switch off GRAPH_PART explicitly ( doing something linke setAndCheckCache )
setAndCheckCache ( METIS GRAPH_PART )
=======
if ( NOT CMAKE_BUILD_TYPE )
   set ( CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif ( NOT CMAKE_BUILD_TYPE )

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
>>>>>>> ac50f33d4c1bacb2d4232fabfbb762a43073ac10:CMake/Modules/SetBuildFlags.cmake

## Check if cache variable is already set
#if    ( DEFINED USE_GRAPH_PART )
#	# do nothing
## if cache variable is NOT set
#else ( DEFINED USE_GRAPH_PART )
#	# Check if package was found
#    if    ( METIS_FOUND ) # ParMetis can only be found with Metis
#    	set ( USE_PACKAGE TRUE )
#    else  ( METIS_FOUND )
#        set ( USE_PACKAGE FALSE )
#    endif ( METIS_FOUND )
#              
#    # Set cache variable
#    set ( ${CACHE_VARIABLE_NAME} ${USE_PACKAGE} CACHE BOOL "Enable / Disable use of ${PACKAGE_NAME}" )
#endif ( DEFINED USE_GRAPH_PART )
