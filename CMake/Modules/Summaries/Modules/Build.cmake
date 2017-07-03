###
 # @file CMake/Modules/Summaries/Modules/Build.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief Summary concerning general build options.
 # @author Lauretta Schubert
 # @date 04.05.2016
###

heading ( "Build options:" "" )

# EXAMPLES
heading3 ( "Examples" "BUILD_EXAMPLES" )

# TEST
heading3 ( "Test" "BOOST_TEST_ENABLED" )
found_message ( "Boost Unit Test" "Boost_UNIT_TEST_FRAMEWORK_FOUND" "OPTIONAL" "Version ${BOOST_VERSION} at ${BOOST_INCLUDE_DIR}" )

# DOC
heading3 ( "Documentation" "DOC_ENABLED" )
    found_message ( "Sphinx" "SPHINX_FOUND" "OPTIONAL" "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )
    found_message ( "Doxygen" "DOXYGEN_FOUND" "OPTIONAL" "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )
