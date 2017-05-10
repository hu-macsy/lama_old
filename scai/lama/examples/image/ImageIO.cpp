/**
 * @file ImageIO.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of routines to read/write image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/GridVector.hpp>

#include <scai/lama/examples/image/ImageIO.hpp>
#include <scai/lama/examples/image/BitmapIO.hpp>
#include <scai/lama/examples/image/PngIO.hpp>

#include<png.h>
#include<fstream>

namespace scai
{

using namespace utilskernel;
using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( ImageIO::logger, "ImageIO" )

/* ------------------------------------------------------------------------------------ */
/*   Read image file, make choice by suffix                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void ImageIO::read( GridVector<ValueType>& imageData, const std::string& inputFileName )
{
    std::string suffix;
    size_t pos = inputFileName.find_last_of( "." );

    if ( pos != std::string::npos )
    {
        suffix = inputFileName.substr( pos );
    }
    else
    {
        COMMON_THROWEXCEPTION( inputFileName << " as image input file has no suffix." )
    }

    LArray<ValueType> data;   // pixel data
    common::Grid3D grid( 0, 0, 0 );  // default grid

    if ( suffix == ".png" )
    {
        PngIO io;
        io.read( data, grid, inputFileName );
        SCAI_LOG_INFO( logger, "read PNG file " << inputFileName << ", data = " << data << ", grid = " << grid )
        imageData.swap( data, grid );
    }
    else if ( suffix == ".bmp" )
    {
        BitmapIO io;
        io.read( data, grid, inputFileName );
        SCAI_LOG_INFO( logger, "read BMP file " << inputFileName << ", data = " << data << ", grid = " << grid )
        imageData.swap( data, grid );
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Write image file, make choice by suffix                                            */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void ImageIO::write( const GridVector<ValueType>& imageData, const std::string& outputFileName )
{
    SCAI_LOG_INFO( logger, "write image data = " << imageData << " to file " << outputFileName )
    std::string suffix;
    size_t pos = outputFileName.find_last_of( "." );

    if ( pos != std::string::npos )
    {
        suffix = outputFileName.substr( pos );
    }
    else
    {
        COMMON_THROWEXCEPTION( outputFileName << " as image output file has no suffix." )
    }

    const HArray<ValueType>& pixelData = imageData.getLocalValues();
    const common::Grid& grid = imageData.globalGrid();
    SCAI_ASSERT_EQ_ERROR( 3, grid.nDims(), "no RGB image data" );
    const common::Grid3D& pixelGrid = reinterpret_cast<const common::Grid3D&>( grid );

    if ( suffix == ".png" )
    {
        PngIO io;
        io.write( pixelData, pixelGrid, outputFileName );
    }
    else if ( suffix == ".bmp" )
    {
        BitmapIO io;
        io.write( pixelData, pixelGrid, outputFileName );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unknown suffix for writing image: " << suffix )
    }
}

// instantiate methods

template
void ImageIO::write( const GridVector<float>&, const std::string& );

template
void ImageIO::read( GridVector<float>&, const std::string& );

} /* end namespace lama */

} /* end namespace scai */
