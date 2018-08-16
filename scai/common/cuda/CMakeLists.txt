###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMakeLists of common/cuda
 # @author Thomas Brandes
 # @date 09.06.2015
###

if ( CUDA_FOUND AND USE_CUDA )
    set ( CUDACLASSES CUDAError CUDASettings CUDACtx CUDAAccess )
    set ( CUDAHEADERS CUDATexVector CUDAUtils )
endif ()

# add source and header files to the current project

scai_project (

    CLASSES

        ${CUDACLASSES}        

    HEADERS

        launchHelper
        CUDACallable
        
        ${CUDAHEADERS}

    ADD_PARENT_SCOPE
)
