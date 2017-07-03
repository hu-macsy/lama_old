###
 # @file scai/logging/test/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
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
 # @brief CMake configuration file for tests of SCAI module logging
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Note: logging has no Boost unit test, but own test script

## make.inc used by Makefile, will be configured

scai_test_scripts( FILES make.inc CONFIGURE )

scai_test_scripts( SCRIPTS       test.sh
                   FILES         Makefile README complexLogging.cpp simpleLogging.cpp
                   CODE_COVERAGE ${USE_CODE_COVERAGE} )

# add test ( cmake commands )

enable_testing()

add_test ( logging_test test.sh )

