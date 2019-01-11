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
 # @brief sparsekernel/cuda/cusparse/CMakeLists.txt
 # @author Eric Schricker
 # @date 18.02.2016
###

if ( CUDA_cusparse_LIBRARY )

### Add classes and header files ###

### Currently CUsparse can only be used if IndexType == int is used

if    ( SCAI_INDEX_TYPE STREQUAL "int" )
    set ( CUSPARSE_CLASSES CUSparseCSRUtils )
    set ( CUSPARSE_HEADERS CUSPARSETrait CUSPARSEWrapper )
endif ( SCAI_INDEX_TYPE STREQUAL "int" )

scai_project ( CUDA_CLASSES ${CUSPARSE_CLASSES} HEADERS ${CUSPARSE_HEADERS} ADD_PARENT_SCOPE )

endif ()
