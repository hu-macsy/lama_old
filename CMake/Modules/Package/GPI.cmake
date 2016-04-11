###
 # @file package/GPI2.cmake
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief findPackage and configuration of GPI
 # @author Thomas Brandes
 # @date 15.02.2016
 # @since 2.0.0
###

### USE_GPI              - if GPI is enabled
### GPI_FOUND            - if GPI2 and IBVERBS are founds
### SCAI_GPI_INCLUDE_DIR - GPI include directory
### SCAI_GPI_LIBRARIES   - all needed GPI libraries (GPI2 and ibverbs)
### GPI_ENABLED          - if GPI_FOUND AND USE_GPI
### GPI_VERSION

find_package ( GPI2 ${SCAI_FIND_PACKAGE_FLAGS} )

if    ( GPI2_FOUND )
    ## get GPI2 version
    try_run ( GPI2_RUN_RESULT_VAR GPI2_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}
        ${CMAKE_MODULE_PATH}/VersionCheck/gpi.cpp
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${GPI2_INCLUDE_DIR}
        COMPILE_OUTPUT_VARIABLE GPI2_COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE GPI2_RUN_OUTPUT_VAR )

    set ( GPI2_VERSION ${GPI2_RUN_OUTPUT_VAR} )
endif ( GPI2_FOUND )

if ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    message ( STATUS "GPI2_FOUND=${GPI2_FOUND}" )
    message ( STATUS "GPI2_INCLUDE_DIR=${GPI2_INCLUDE_DIR}" )
    message ( STATUS "GPI2_LIBRARIES=${GPI2_LIBRARIES}" )
endif ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )

find_package ( Ibverbs ${SCAI_FIND_PACKAGE_FLAGS} )

if ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    message ( STATUS "IBVERBS_FOUND=${IBVERBS_FOUND}" )
    message ( STATUS "IBVERBS_INCLUDE_DIR=${IBVERBS_INCLUDE_DIR}" )
    message ( STATUS "IBVERBS_LIBRARIES=${IBVERBS_LIBRARIES}" )
endif ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )

### ALLOW to switch off GPI2 explicitly ###
# do what setAndCheckCache does but with 2 packages
# Check if cache variable is already set
if    ( DEFINED USE_GPI )
    # do nothing
    # if cache variable is NOT set
    set( USE_GPI ${USE_GPI} )

else ( DEFINED USE_GPI )
    # Check if package was found
    if    ( GPI2_FOUND AND IBVERBS_FOUND )
        set ( USE_PACKAGE TRUE )
    else  ( GPI2_FOUND AND IBVERBS_FOUND )
        set ( USE_PACKAGE FALSE )
    endif ( GPI2_FOUND AND IBVERBS_FOUND )
        
    # Set cache variable
    set ( USE_GPI ${USE_PACKAGE} )

endif ( DEFINED USE_GPI )

set ( GPI_FOUND FALSE )
set ( GPI_ENABLED FALSE )

if    ( USE_GPI AND GPI2_FOUND AND IBVERBS_FOUND )

    # message ( STATUS "set SCAI_GPI_INCLUDE_DIR + SCAI_GPI_LIBRARIES" )

    # just for making it the same variable ending for all packages
    set ( SCAI_GPI_INCLUDE_DIR ${GPI2_INCLUDE_DIR} ${IBVERBS_INCLUDE_DIR} )
    
    # conclude all needed CUDA libraries
    set ( SCAI_GPI_LIBRARIES ${GPI2_LIBRARIES} ${IBVERBS_LIBRARIES} )

    set ( GPI_FOUND TRUE )
    set ( GPI_ENABLED TRUE )

endif ( USE_GPI AND GPI2_FOUND AND IBVERBS_FOUND )

