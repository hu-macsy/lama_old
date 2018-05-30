###
 # @file CMakeLists.txt
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
 # @brief lama/storage/CMakeLists.txt
 # @author Jan Ecker
 # @date 28.03.2013
###

### Add classes and header files ###

set ( CLASSES

        Format

        _MatrixStorage
        MatrixStorage

        COOStorage
        CSRStorage
        DenseStorage
        DIAStorage
        ELLStorage
        JDSStorage
        StencilStorage

        StorageMethods
    )
    
scai_project ( CLASSES ${CLASSES} HEADERS ${HEADERS} ADD_PARENT_SCOPE )
