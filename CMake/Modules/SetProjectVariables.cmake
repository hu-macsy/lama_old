###
 # @file Variables.cmake
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
 # @brief set project name variable for flexible use in the subprojects
 # @author Lauretta Schubert
 # @date 04.02.2016
 # @since 2.0.0
###

## upper case project name
string ( TOUPPER ${PROJECT_NAME} UPPER_PROJECT_NAME )

## get project "sur"name (everthing after scai_)
string ( SUBSTRING ${PROJECT_NAME} 5 -1 PROJECT_SURNAME )

## upper case project name for project dependent variables
set ( ${UPPER_PROJECT_NAME}_INCLUDE_DIR include/scai/${PROJECT_SURNAME} )

## example dir
set ( PROJECT_EXAMPLE_DIR "share/examples/scai-${PROJECT_SURNAME}-${${UPPER_PROJECT_NAME}_VERSION}" )

## doc dir
set ( DOC_ROOT_DIR share/doc )
set ( SPHINX_ROOT_DIR "${DOC_ROOT_DIR}/scai-${SCAI_LAMA_ALL_VERSION}")
set ( PROJECT_DOC_DIR "${SPHINX_ROOT_DIR}/scai-${PROJECT_SURNAME}-${${UPPER_PROJECT_NAME}_VERSION}" )
