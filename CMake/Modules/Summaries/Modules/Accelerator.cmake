###
 # @file CMake/Modules/Summaries/Modules/Accelerator.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief Summary concerning Accelerator support.
 # @author Lauretta Schubert
 # @date 08.04.2016
###

# OpenMP usage
heading3 ( "OpenMP" "USE_OPENMP" )
    found_message ( "OpenMP" "OPENMP_VERSION" "OPTIONAL" "Version ${OPENMP_VERSION}" )
    found_message ( "compile flag" "OpenMP_CXX_FLAGS" "OPTIONAL" "${OpenMP_CXX_FLAGS}" )
    found_message ( "schedule type" "SCAI_OMP_SCHEDULE" "OPTIONAL" "set to \"${SCAI_OMP_SCHEDULE}\"" )

# LAMA CUDA
heading3 ( "CUDA" "CUDA_ENABLED" )
    found_message ( "CUDA" "CUDA_FOUND" "OPTIONAL" "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
    found_message ( "Compute Capability" "CUDA_HAVE_GPU" "OPTIONAL" "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
heading3 ( "MIC" "USE_MIC" )