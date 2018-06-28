###
 # @file blaskernel/cuda/CMakeLists.txt
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
 # @brief CMake configuration for CUDA implementations of BLAS kernels.
 # @author Eric Schricker
 # @date 11.11.2015
###

### Add classes and header files ###

scai_project(

    CLASSES      # .cpp, .hpp

        CUDABLAS2
        CUDABLAS3

    CUDA_CLASSES  # .cu, .hpp

        CUDABLAS1

    CUDA_SOURCES   # .cu


    HEADERS        # .hpp

        CUBLASWrapper
	    CUBLASTrait

    ADD_PARENT_SCOPE
)
