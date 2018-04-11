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

#include <scai/lama/io/PngIO.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/Grid.hpp>

#include<png.h>
#include<fstream>

#define MAT_SUFFIX ".png" 

namespace scai
{

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( PngIO::logger, "FileIO.Png" )

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* PngIO::create()
{
    return new PngIO();
}

std::string PngIO::createValue()
{
    return MAT_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool PngIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* ------------------------------------------------------------------------------------ */

PngIO::PngIO()
{
    SCAI_LOG_INFO( logger, "PngIO uses libpng version " << png_libpng_ver );
}

/* ------------------------------------------------------------------------------------ */
/*   Read PNG File                                                                   */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void PngIO::readGridImpl( HArray<ValueType>& imageData, common::Grid& imageSize, const std::string& inputFileName )
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


    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
        COMMON_THROWEXCEPTION( "Serious error during reading png file " << inputFileName )
    }

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
void PngIO::writeGridImpl( const HArray<ValueType>& data, const common::Grid& grid, const std::string& outputFileName )
{
    SCAI_LOG_INFO( logger, "write image, shape = " << grid )
    const IndexType height = grid.size( 0 );
    const IndexType width  = grid.size( 1 );
    const IndexType ncolor = grid.size( 2 );
    SCAI_ASSERT_EQ_ERROR( 3, ncolor, "3 values per pixel expected" )
    // convert the array data to png data
    png_bytep* row_pointers = new png_bytep[ static_cast<int>( height ) ];
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


    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        COMMON_THROWEXCEPTION( "Serious error during writing png file " << outputFileName )
    }


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

void PngIO::readGridArray( _HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    IOWrapper<PngIO, SCAI_TYPELIST( float, double )>::readGridImpl( ( PngIO& ) *this, data, grid, inputFileName );
}

void PngIO::writeGridArray( const _HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    IOWrapper<PngIO, SCAI_TYPELIST( float, double )>::writeGridImpl( ( PngIO& ) *this, data, grid, outputFileName );
}

/* --------------------------------------------------------------------------------- */

void PngIO::writeStorage( const _MatrixStorage& storage, const std::string& outputFileName )
{
    // todo: write dense storage as bitmap, maybe greyscale 

    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "write storage " << storage << " to " << outputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::readStorage(
    _MatrixStorage& storage,
    const std::string& inputFileName,
    const IndexType offsetRow,
    const IndexType nRows )
{
    storage.clear();

    SCAI_ASSERT_EQ_ERROR( 0, offsetRow, "No chunk read for bitmap file" )
    SCAI_ASSERT_EQ_ERROR( invalidIndex, nRows, "No chunk read for bitmap file" )

    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: read storage from " << inputFileName )
}

/* --------------------------------------------------------------------------------- */

std::string PngIO::getMatrixFileSuffix() const
{
    return createValue();
}

/* --------------------------------------------------------------------------------- */

std::string PngIO::getVectorFileSuffix() const
{
    return createValue();
}

/* --------------------------------------------------------------------------------- */

void PngIO::writeArray( const hmemo::_HArray& array, const std::string& outputFileName )
{
    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: write array " << array << " to " << outputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::writeSparse( 
    const IndexType n, 
    const hmemo::HArray<IndexType>& indexes, 
    const hmemo::_HArray& values, 
    const std::string& outputFileName )
{
    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: write spare array ( n = " << n 
                         << ", indexes = " << indexes << ", values = " << values << " ) to " << outputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& inputFileName )
{
    numRows    = 0;
    numColumns = 0;
    numValues  = 0;

    COMMON_THROWEXCEPTION( "Unsupported for bitmap file: read storage info from file " << inputFileName )
}

void PngIO::readArrayInfo( IndexType& size, const std::string& inputFileName )
{
    size = 0;

    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: read array info from file " << inputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::readArray( hmemo::_HArray& array, const std::string& inputFileName, const IndexType offset, const IndexType n )
{
    array.clear();

    SCAI_ASSERT_EQ_ERROR( 0, offset, "chunk read not supported" )
    SCAI_ASSERT_EQ_ERROR( n, invalidIndex, "chunk read not supported" )

    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: read array ( offset = " << offset << ", n = " << " ) from file " << inputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const std::string& inputFileName )
{
    size = 0;
    indexes.clear();
    values.clear();

    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: read sparse array from " << inputFileName )
}

/* --------------------------------------------------------------------------------- */

void PngIO::writeAt( std::ostream& stream ) const
{
    stream << "PngIO ( suffix = " << MAT_SUFFIX << ", ";
    stream << ", only image data )";
}

/* --------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
