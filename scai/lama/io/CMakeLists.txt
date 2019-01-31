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
 # @brief lama/io/CMakeLists.txt
 # @author Jan Ecker
 # @date 02.04.2013
###

### Add classes and header files ###

set ( IO_CLASSES 

        IOStream

        FileIO
        SAMGIO
        LamaIO
        MatrixMarketIO
        PETScIO
        TextIO
        PartitionIO
        ImageIO
        BitmapIO
) 

## IO for Matlab format only possible if ZLIB is available

if ( ZLIB_FOUND AND USE_ZLIB )
    set ( IO_CLASSES ${IO_CLASSES} MatlabIO MATIOStream )
endif ()

## IO for PNG format only possible if libpng is available 

if ( PNG_FOUND AND USE_PNG )
    set ( IO_CLASSES ${IO_CLASSES} PngIO )
endif ()

scai_project (

    CLASSES
        
        ${IO_CLASSES}

    HEADERS

        IOWrapper

    ADD_PARENT_SCOPE
)
