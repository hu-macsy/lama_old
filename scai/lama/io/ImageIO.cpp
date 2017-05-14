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
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/lama/io/ImageIO.hpp>
#include <scai/lama/io/BitmapIO.hpp>
#include <scai/lama/io/PngIO.hpp>

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
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported suffix (" << suffix << ") for reading image file, only .bmp or .png" )
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
        COMMON_THROWEXCEPTION( "Unsupported suffix (" << suffix << ") for writing image file, only .bmp or .png" )
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Write scaled image                                                                 */
/* ------------------------------------------------------------------------------------ */

static double colorMap[] =
{
   0.00000,  0.00000,  0.56250,
   0.00000,  0.00000,  0.62500,
   0.00000,  0.00000,  0.68750,
   0.00000,  0.00000,  0.75000,
   0.00000,  0.00000,  0.81250,
   0.00000,  0.00000,  0.87500,
   0.00000,  0.00000,  0.93750,
   0.00000,  0.00000,  1.00000,
   0.00000,  0.06250,  1.00000,
   0.00000,  0.12500,  1.00000,
   0.00000,  0.18750,  1.00000,
   0.00000,  0.25000,  1.00000,
   0.00000,  0.31250,  1.00000,
   0.00000,  0.37500,  1.00000,
   0.00000,  0.43750,  1.00000,
   0.00000,  0.50000,  1.00000,
   0.00000,  0.56250,  1.00000,
   0.00000,  0.62500,  1.00000,
   0.00000,  0.68750,  1.00000,
   0.00000,  0.75000,  1.00000,
   0.00000,  0.81250,  1.00000,
   0.00000,  0.87500,  1.00000,
   0.00000,  0.93750,  1.00000,
   0.00000,  1.00000,  1.00000,
   0.06250,  1.00000,  0.93750,
   0.12500,  1.00000,  0.87500,
   0.18750,  1.00000,  0.81250,
   0.25000,  1.00000,  0.75000,
   0.31250,  1.00000,  0.68750,
   0.37500,  1.00000,  0.62500,
   0.43750,  1.00000,  0.56250,
   0.50000,  1.00000,  0.50000,
   0.56250,  1.00000,  0.43750,
   0.62500,  1.00000,  0.37500,
   0.68750,  1.00000,  0.31250,
   0.75000,  1.00000,  0.25000,
   0.81250,  1.00000,  0.18750,
   0.87500,  1.00000,  0.12500,
   0.93750,  1.00000,  0.06250,
   1.00000,  1.00000,  0.00000,
   1.00000,  0.93750,  0.00000,
   1.00000,  0.87500,  0.00000,
   1.00000,  0.81250,  0.00000,
   1.00000,  0.75000,  0.00000,
   1.00000,  0.68750,  0.00000,
   1.00000,  0.62500,  0.00000,
   1.00000,  0.56250,  0.00000,
   1.00000,  0.50000,  0.00000,
   1.00000,  0.43750,  0.00000,
   1.00000,  0.37500,  0.00000,
   1.00000,  0.31250,  0.00000,
   1.00000,  0.25000,  0.00000,
   1.00000,  0.18750,  0.00000,
   1.00000,  0.12500,  0.00000,
   1.00000,  0.06250,  0.00000,
   1.00000,  0.00000,  0.00000,
   0.93750,  0.00000,  0.00000,
   0.87500,  0.00000,  0.00000,
   0.81250,  0.00000,  0.00000,
   0.75000,  0.00000,  0.00000,
   0.68750,  0.00000,  0.00000,
   0.62500,  0.00000,  0.00000,
   0.56250,  0.00000,  0.00000,
   0.50000,  0.00000,  0.00000
};

/** Interpolate a color value between two colors.
 *
 *  @param[in]  color1 array with rgb values of 1st color
 *  @param[in]  color2 array with rgb values of 2nd color
 *  @param[in]  factor 0.0 stands for color1, 1.0 for color2
 *  @param[out] color array will contain the interpolated rgb values
 */
static void interpolateColor( double color[3], const double color1[3], const double color2[3], const double correction )
{
    color[0] = ( 1.0 - correction ) * color1[0] + correction * color2[0];
    color[1] = ( 1.0 - correction ) * color1[1] + correction * color2[1];
    color[2] = ( 1.0 - correction ) * color1[2] + correction * color2[2];
}

/** Compute for a value in range [0:1] a corresponding color entry by a color map 
 *
 *  @param[in] value must be between 0.0 and 1.0
 *  @param[out] color is the computed color 
 */

static void computeColor( double color[3], const double value )
{
    static IndexType nColorsInTable = sizeof( colorMap ) / sizeof( double ) / 3;

    static IndexType nSections = nColorsInTable - 1;
    static double    width     = double( 1 ) / double( nSections );

    // make sure that value is in range [0,1]

    double v = common::Math::min( 1.0, value );
    v = common::Math::max( 0.0, v );

    // find entry 

    IndexType colorIndex = common::Math::floor( value * nSections );

    if ( colorIndex >= nSections )
    {
        colorIndex = nSections - 1;
    }


    // now scale color betwen two neighbored colors

    double correction = ( value - double( colorIndex ) / nSections ) / width;

    interpolateColor( color, &colorMap[ 3 * colorIndex], &colorMap[ 3 * colorIndex + 3], correction );
}

template<typename ValueType>
void ImageIO::writeSC( const GridVector<ValueType>& arrayData, const std::string& outputFileName )
{
    SCAI_LOG_ERROR( logger, "write scaled array data = " << arrayData << " to file " << outputFileName )

    SCAI_ASSERT_EQ_ERROR( 2, arrayData.nDims(), "writeSC only for two-dimensional arrays" )

    const common::Grid2D& arrayGrid = reinterpret_cast<const common::Grid2D&>( arrayData.globalGrid() );

    common::Grid3D imageGrid( arrayGrid.size(0), arrayGrid.size(1), 3 );

    GridVector<float> imageData( imageGrid, 0 );

    Scalar min = arrayData.min();
    SCAI_LOG_ERROR( logger, "array min = " << min )
    ValueType minVal = min.getValue<ValueType>();
    Scalar max = arrayData.max();
    SCAI_LOG_ERROR( logger, "array max = " << min )
    ValueType maxVal = max.getValue<ValueType>();

    if ( minVal == maxVal )
    {
        minVal -= 0.01;
        maxVal += 0.01;
    }

    ValueType scale = maxVal - minVal;

    const IndexType nColorsInTable = sizeof( colorMap ) / sizeof( double ) / 3;

    SCAI_LOG_ERROR( logger, "Colormap has " << nColorsInTable << " entries, scale range = " << minVal << " - " << maxVal )

    {
        GridReadAccess<ValueType> rArray( arrayData );
        GridWriteAccess<float> wImage( imageData );

        for ( IndexType i = 0; i < arrayGrid.size( 0 ); ++i )
        {
            for ( IndexType j = 0; j < arrayGrid.size( 1 ); ++j )
            {
                ValueType val = rArray( i, j );
                ValueType scaledVal = ( val - minVal ) / scale;  // val in [0,1]

                double color[3];

                computeColor( color, scaledVal ); 

                wImage( i, j, 0 ) = color[ 0 ] * 255.0 + 0.5;
                wImage( i, j, 1 ) = color[ 1 ] * 255.0 + 0.5;
                wImage( i, j, 2 ) = color[ 2 ] * 255.0 + 0.5;
            }
        }
    }

    ImageIO::write( imageData, outputFileName );
}

// instantiate methods for supported array/vector types

#define SCAI_IMAGE_IO_INSTANTIATIONS( _type )                              \
                                                                           \
    template COMMON_DLL_IMPORTEXPORT                                       \
    void ImageIO::read( GridVector<_type>&, const std::string& );          \
                                                                           \
    template COMMON_DLL_IMPORTEXPORT                                       \
    void ImageIO::write( const GridVector<_type>&, const std::string& );   \
                                                                           \
    template COMMON_DLL_IMPORTEXPORT                                       \
    void ImageIO::writeSC( const GridVector<_type>&, const std::string& ); \

SCAI_COMMON_LOOP( SCAI_IMAGE_IO_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
