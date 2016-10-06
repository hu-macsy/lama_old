###
 # @file CMake/Modules/Settings/switchChoices.cmake
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
 # @brief Sets all lists of choices that can be switch in the LAMA configutation and defines their default behaviour.
 # @author Lauretta Schubert
 # @date 01.04.2016
###

set ( TRUE_FALSE_CHOICES ON OFF )

# Possible and default choices for numeric types and IndexType 

set ( NONE_COMPLEX_VALUES "float" "double" "long double" )
set ( COMPLEX_VALUES      "ComplexFloat" "ComplexDouble" "ComplexLongDouble" )
set ( TYPE_LIST_VALUES    ${NONE_COMPLEX_VALUES} ${COMPLEX_VALUES} )
set ( LONG_DOUBLE_VALUES  "long double" "ComplexLongDouble" )
set ( INDEX_TYPE_OPTIONS  "int" "long" "unsigned int" "unsigned long" )
set ( INDEX_TYPE_DEFAULT  "int" )

set ( BUILD_DOC_DEFAULT ON )
set ( BUILD_EXAMPLES_DEFAULT ON )
set ( BUILD_TEST_DEFAULT ON )
set ( USE_CODE_COVERAGE_DEFAULT OFF )

set ( CMAKE_BUILD_TYPE_CHOICES "None" "Debug" "Release" "RelWithDebInfo" "MinSizeRel" ) 
set ( CMAKE_BUILD_TYPE_DEFAULT "Debug" )

set ( SCAI_ASSERT_CHOICES DEBUG ERROR OFF )
set ( SCAI_ASSERT_DEFAULT DEBUG )

set ( SCAI_BLAS_LIBRARY_CHOICES auto MKL BLAS INTERNALBLAS )
set ( SCAI_BLAS_LIBRARY_DEFAULT auto )

set ( SCAI_DOC_TYPE_CHOICES html json xml latex )
set ( SCAI_DOC_TYPE_DEFAULT html )

set ( SCAI_LIBRARY_TYPE_CHOICES STATIC SHARED )
set ( SCAI_LIBRARY_TYPE_DEFAULT SHARED )

set ( SCAI_LOGGING_CHOICES TRACE DEBUG INFO WARN ERROR OFF )
# no default, decision depending on choosen CMAKE_BUILD_TYPE --> see Settings/logging
