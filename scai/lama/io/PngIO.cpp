/**
 * @file PngIO.cpp
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

#include <scai/lama/io/PngIO.hpp>

#include<png.h>
#include<fstream>

namespace scai
{

using namespace utilskernel;
using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( PngIO::logger, "ImageIO.Png" )

/* ------------------------------------------------------------------------------------ */
/*   Read PNG File                                                                   */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void PngIO::readImpl( HArray<ValueType>& imageData, common::Grid& imageSize, const std::string& inputFileName )
{
    char header[8];    // 8 is the maximum size that can be checked
    /* open file and test for it being a png */
    FILE* fp = fopen( inputFileName.c_str(), "rb" );
    SCAI_ASSERT_ERROR( fp, inputFileName << " could not be opened for reading " )
    size_t readBytes = fread( header, 1, 8, fp );
    SCAI_ASSERT_EQ_ERROR( 8, readBytes, "Could not read header" )

    if ( png_sig_cmp( reinterpret_cast<png_byte*>( header ), 0, 8 ) )
    {
        COMMON_THROWEXCEPTION( inputFileName << " is not recognized as a PNG file" )
    }

    /* read basic infomation */

    png_structp png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
    SCAI_ASSERT_ERROR( png_ptr, "png_create_read_struct failed" )

    png_infop info_ptr = png_create_info_struct( png_ptr );
    SCAI_ASSERT_ERROR( info_ptr, "png_create_info_struct failed" )

    // prepare error handling, i.e. define code that is used for error handling

    png_init_io( png_ptr, fp );
    png_set_sig_bytes( png_ptr, 8 );
    png_read_info( png_ptr, info_ptr );
    IndexType width = png_get_image_width( png_ptr, info_ptr );
    IndexType height = png_get_image_height( png_ptr, info_ptr );

    SCAI_LOG_INFO( logger, "read image of " << height << " x " << width << " ( height x width )" )

    png_byte color_type = png_get_color_type( png_ptr, info_ptr );
    png_byte png_bit_depth = png_get_bit_depth( png_ptr, info_ptr );

    int number_of_passes = png_set_interlace_handling( png_ptr );

    SCAI_LOG_INFO( logger, "color_type = " << static_cast<int>( color_type )
                    << ", png_bit_depth = " << static_cast<int>( png_bit_depth )
                    << ", passes = " << number_of_passes )

    SCAI_ASSERT_EQ_ERROR( 8, png_bit_depth, "read only 8-bit depth supported" )

    if ( number_of_passes != 1 )
    {
        SCAI_LOG_WARN( logger, "Reading " << inputFileName << ": multiple passes not supported"
                               << ", #passes = " << number_of_passes )
    }

    png_read_update_info( png_ptr, info_ptr );

    /* read file */

    png_bytep* row_pointers = new png_bytep[ height ];

    for ( IndexType y = 0; y < height; y++ )
    {
        row_pointers[y] = new png_byte[ png_get_rowbytes( png_ptr, info_ptr ) ];
    }

    png_read_image( png_ptr, row_pointers );

    fclose( fp );

    IndexType valsPerPixel;

    switch ( color_type )
    {  
        case PNG_COLOR_TYPE_RGBA :  valsPerPixel = 4; 
                                    break;
        case PNG_COLOR_TYPE_RGB  :  valsPerPixel = 3;
                                    break;

        default: COMMON_THROWEXCEPTION( "color type = " << ( int ) color_type << " not supported" )
    }

    imageSize = common::Grid3D( height, width, 3 );
    {
        WriteOnlyAccess<ValueType> wX( imageData, imageSize.size() );

        for ( IndexType i = 0; i < height; ++i )
        {
            png_bytep row_ptr = row_pointers[ i ];

            for ( IndexType j = 0; j < width; ++j )
            {
                wX[ 3 * width * i + 3 * j ] = row_ptr[ valsPerPixel * j ];
                wX[ 3 * width * i + 3 * j + 1 ] = row_ptr[ valsPerPixel * j + 1 ];
                wX[ 3 * width * i + 3 * j + 2 ] = row_ptr[ valsPerPixel * j + 2 ];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Write Png file                                                                     */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void PngIO::writeImpl( const HArray<ValueType>& data, const common::Grid& grid, const std::string& outputFileName )
{
    SCAI_LOG_INFO( logger, "write image, shape = " << grid )
    const IndexType height = grid.size( 0 );
    const IndexType width  = grid.size( 1 );
    const IndexType ncolor = grid.size( 2 );
    SCAI_ASSERT_EQ_ERROR( 3, ncolor, "3 values per pixel expected" )
    // convert the array data to png data
    png_bytep* row_pointers = new png_bytep[ height ];
    {
        ReadAccess<ValueType> rX( data );

        for ( IndexType y = 0; y < height; y++ )
        {
            row_pointers[y] = new png_byte[ 4 * width ];

            for ( IndexType x = 0; x < width; ++x )
            {
                // png_bitmap requires r, g, b, opaque value
                row_pointers[y][ 4 * x ]     = static_cast<png_byte>( rX[ 3 * width * y + 3 * x ] );
                row_pointers[y][ 4 * x + 1 ] = static_cast<png_byte>( rX[ 3 * width * y + 3 * x + 1 ] );
                row_pointers[y][ 4 * x + 2 ] = static_cast<png_byte>( rX[ 3 * width * y + 3 * x + 2 ] );
                row_pointers[y][ 4 * x + 3 ] = static_cast<png_byte>( 255 );
            }
        }
    }
    png_structp png_ptr;
    FILE* fp = fopen( outputFileName.c_str(), "wb" );
    SCAI_ASSERT_ERROR( fp, "File " << outputFileName << " could not be opened for writing" )
    /* initialize stuff */
    png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
    SCAI_ASSERT_ERROR( png_ptr, "png_create_write_struct failed" );
    png_infop info_ptr = png_create_info_struct( png_ptr );
    SCAI_ASSERT_ERROR( info_ptr, "png_create_info_struct failed" );

    // prepare error handling, i.e. define code that is used for error handling


    png_init_io( png_ptr, fp );

    /* write header */
    png_byte bit_depth = 8;
    png_byte color_type = PNG_COLOR_TYPE_RGBA;
    png_set_IHDR( png_ptr, info_ptr, width, height,
                  bit_depth, color_type, PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );
    png_write_info( png_ptr, info_ptr );

    png_write_image( png_ptr, row_pointers );

    png_write_end( png_ptr, NULL );

    fclose( fp );
}

/* ------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------ */

void PngIO::read( _HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    ImageIOWrapper<PngIO, SCAI_TYPELIST( float, double )>::read( ( PngIO& ) *this, data, grid, inputFileName );
}

void PngIO::write( const _HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    ImageIOWrapper<PngIO, SCAI_TYPELIST( float, double )>::write( ( PngIO& ) *this, data, grid, outputFileName );
}

/* ------------------------------------------------------------------------------------ */

PngIO::PngIO()
{
    SCAI_LOG_INFO( logger, "PngIO uses libpng version " << png_libpng_ver );
}


} /* end namespace lama */

} /* end namespace scai */
