###
 # @file Sphinx.cmake
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
 # @brief find package Sphinx and set BUILD_DOC flag
 # @author Jan Ecker
 # @date 11.11.2015
 # @since 1.0.1
###

find_package ( Sphinx ${SCAI_FIND_PACKAGE_FLAGS} )

# Check if cache variable is already set
if ( DEFINED BUILD_DOC )
    # if use of package is enabled
    if ( ${BUILD_DOC} )
        if ( ${SPHINX_FOUND} )
        else ( ${SPHINX_FOUND} )
            # if package is enabled, but not found: ERROR!
            message ( STATUS "Sphinx missing, but build of doc is enabled!" )
        endif ( ${SPHINX_FOUND} )
    endif ( ${BUILD_DOC} )
# if cache variable is NOT set
else ( DEFINED BUILD_DOC )
    # Check if package was found
    if ( ${SPHINX_FOUND} )
        set ( USE_PACKAGE TRUE )
    else ( ${SPHINX_FOUND} )
        set ( USE_PACKAGE FALSE )
    endif ( ${SPHINX_FOUND} )
    
    # Set cache variable
    set ( BUILD_DOC ${USE_PACKAGE} CACHE BOOL "Enable / Disable building of doc" )
endif ( DEFINED BUILD_DOC )