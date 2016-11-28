###
 # @file Dependencies/minimalSupportedVersions.cmake
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
 # @brief Defines minimal supported library versions for compatiblity check at find_package
 # @author Lauretta Schubert
 # @date 27.10.2016
###

#set ( BOOST_MINIMUM_VERSION 1.36 )
set ( BOOST_TEST_MINIMUM_VERSION 1.41 )

set ( OMP_MINIMUM_VERSION 3.0 ) # because of use of collapse

set ( CUDA_MINIMUM_VERSION 4.0 ) # because we explicitly use cublas_v2
# cusparse since 3.2, cusolver since 7.0
