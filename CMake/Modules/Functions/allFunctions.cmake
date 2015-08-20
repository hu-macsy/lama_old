###
 # @file Functions.cmake
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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

# defined functions:
#     setAndCheckCache: Function for setting USE_{PACKAGE_NAME} variables depending on {PACKAGE_NAME}_FOUND.
include ( Functions/setAndCheckCache )
#     get_relative_path: returns the relative path to the actual directory to the CMAKE_SOURCE_DIR
include ( Functions/getRelativePath )
#     lama_get_relative_path: returns the relative path to the actual directory to the LAMA_SOURCE_DIR (Path of the actual target)
include ( Functions/lamaGetRelativePath )
#     scai_status_message: prints colored text messages
include ( Functions/scaiStatusMessage )
#     checkValue: checks whether the given value is in the value list ( pass list as "${LIST}" (doublequotes !!!) )
include ( Functions/checkValue )
#     scai_generate_blanks: Simple internal helper function that generates a blank string that fits the size of an given STRING to LENGTH
include ( Functions/scaiGenerateBlanks )
#     check_whitelist: checks white list for containing entry and add entry to argument list ( pass whitelist and arglist as "${LIST}" (doublequotes !!!) )
include ( Functions/checkWhiteList )

# defined makros:
#     lama_set_source_dir: sets the LAMA_SOURCE_DIR (used to mark the path of the actual build target
#     lama_classes: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures the installation of the header files
#     lama_sources: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures the installation of the header files
#     lama_headers: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures # the installation of the header files
#     lama_add: Publishes sources and headers in the parent scope
include ( Functions/lamaSourceSolution )
#     scai_summary_message: generates messages for lama summary page
include ( Functions/scaiSummaryMessage )
#     list_contains: checks if value is part of a list
include ( Functions/listContains )
