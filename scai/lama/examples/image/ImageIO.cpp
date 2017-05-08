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
/*   Data structures for BITMAP file                                                    */
/* ------------------------------------------------------------------------------------ */

#pragma pack(2)

struct BITMAPFILEHEADER
{
    unsigned short int bfType;
    unsigned int bfSize;
    unsigned int bfReserved;
    unsigned int bfOffBits;
};

struct BITMAPINFOHEADER
{
    unsigned int biSize;
    int biWidth;
    int biHeight;
    unsigned short int biPlanes;
    unsigned short int biBitCount;
    unsigned int biCompression;
    unsigned int biSizeImage;
    int biXPelsPerMeter;
    int biYPelsPerMeter;
    unsigned int biClrUsed;
    unsigned int biClrImportant;
};

#pragma pack()

struct COLORMASKS
{
    unsigned int red;
    unsigned int green;
    unsigned int blue;
};

struct COLORTABLE
{
    unsigned int* colors;
};

/* ------------------------------------------------------------------------------------ */
/*   Read PNG file                                                                      */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void ImageIO::readBMP( HArray<ValueType>& data, common::Grid3D& grid, const std::string& inputFileName )
{
    FILE* file;
    file = fopen( inputFileName.c_str(), "rb" );
    SCAI_ASSERT( file, "Could not open file " << inputFileName )
    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;
    // COLORMASKS colorMasks;
    // COLORTABLE colorTable;
    size_t nRecords = fread( &fileHeader, sizeof( fileHeader ), 1, file );

    SCAI_ASSERT_EQ_ERROR( 1, nRecords, "failed to read header" )
    SCAI_ASSERT_EQ_ERROR( fileHeader.bfType, 0x4D42, "not bitmap file" )


    SCAI_LOG_INFO( logger, "fileHeader, bfSize = " << fileHeader.bfSize
                           << ", bfReserved = " << fileHeader.bfReserved 
                           << ", bfOffBits = " << fileHeader.bfOffBits )


    nRecords = fread( &infoHeader, sizeof( infoHeader ), 1, file );

    SCAI_LOG_INFO( logger, "infoHeader, biSize = " << infoHeader.biSize )
    SCAI_LOG_INFO( logger, "infoHeader, biCompression = " << infoHeader.biCompression )
    SCAI_LOG_INFO( logger, "infoHeader, biPlanes = " << infoHeader.biPlanes )
    SCAI_LOG_INFO( logger, "infoHeader, biBitCount = " << infoHeader.biBitCount )
    SCAI_LOG_INFO( logger, "infoHeader, biCompression = " << infoHeader.biCompression )
    SCAI_LOG_INFO( logger, "infoHeader, biSizeImage = " << infoHeader.biSizeImage )
    SCAI_LOG_INFO( logger, "infoHeader, biXPelsPerMeter = " << infoHeader.biXPelsPerMeter )
    SCAI_LOG_INFO( logger, "infoHeader, biYPelsPerMeter = " << infoHeader.biYPelsPerMeter )
    SCAI_LOG_INFO( logger, "infoHeader, biClrUsed = " << infoHeader.biClrUsed )
    SCAI_LOG_INFO( logger, "infoHeader, biClrImportant = " << infoHeader.biClrImportant )

    IndexType imageSize = infoHeader.biSizeImage;

    IndexType width = infoHeader.biWidth;
    IndexType height = infoHeader.biHeight;
    SCAI_LOG_INFO( logger, "BMP file, imag is " << width << " x " << height )

    common::scoped_array<unsigned char> tmpData( new unsigned char[ 3 * width * height] );

    if ( infoHeader.biCompression == 3 )
    {
        // file.read( colorMasks, sizeof( colorMasks ) );
    }

    if ( infoHeader.biBitCount == 24 )
    {
        grid = common::Grid3D( height, width, 3 );
        SCAI_ASSERT_EQ_ERROR( imageSize, grid.size(), "size mismatch" )
        common::scoped_array<unsigned char> tmp( new unsigned char[ grid.size() ] );
        IndexType nBytes = fread( tmp.get(), 1, grid.size(), file );
        SCAI_ASSERT_EQ_ERROR( grid.size(), nBytes, "could not read all pixels" )
        WriteOnlyAccess<ValueType> wData( data, grid.size() );
        {
            for ( IndexType i = 0; i < height; ++i )
            {
                for ( IndexType j = 0; j < width; ++j )
                {
                    // BMP has data from top to bottom, but we need the inverse
                    IndexType bmpPos = 3 * ( i * width + j );
                    IndexType imgPos = 3 * ( ( height - 1 - i ) * width + j );
                    wData[ imgPos   ] = static_cast<ValueType>( tmp[bmpPos + 2] );
                    wData[ imgPos + 1 ] = static_cast<ValueType>( tmp[bmpPos + 1] );
                    wData[ imgPos + 2 ] = static_cast<ValueType>( tmp[bmpPos] );
                }
            }
        }
    }
    else if ( infoHeader.biBitCount == 32 )
    {
        COMMON_THROWEXCEPTION( "32-bit not supported" )
    }
}

template<typename ValueType>
void ImageIO::writeBMP( const HArray<ValueType>& data, const common::Grid3D& grid, const std::string& outputFileName )
{
    FILE* f = fopen( outputFileName.c_str(), "wb" );

    IndexType height  = grid.size( 0 );
    IndexType width   = grid.size( 1 );

    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;

    unsigned int headerBytes = sizeof( fileHeader ) + sizeof( infoHeader );

    fileHeader.bfType     = 0x4D42;
    fileHeader.bfSize     = grid.size() + headerBytes;
    fileHeader.bfReserved = 0;
    fileHeader.bfOffBits  = headerBytes;

    infoHeader.biSize = sizeof( infoHeader );
    infoHeader.biWidth = width;
    infoHeader.biHeight = height;
    infoHeader.biPlanes = 1;
    infoHeader.biBitCount = 24;
    infoHeader.biCompression = 0;
    infoHeader.biSizeImage = grid.size();
    infoHeader.biXPelsPerMeter = 13780;
    infoHeader.biYPelsPerMeter = 13780;
    infoHeader.biClrUsed = 0;
    infoHeader.biClrImportant = 0;

    fwrite( &fileHeader, 1, sizeof( fileHeader ), f );
    fwrite( &infoHeader, 1, sizeof( infoHeader ), f );

    common::scoped_array<unsigned char> tmpData( new unsigned char[ 3 * width * height ] );

    ReadAccess<ValueType> rPixelData( data );

    for ( IndexType i = 0; i < height; ++i )
    {
        for ( IndexType j = 0; j < width; ++j )
        {
            // BMP has data from top to bottom, but we need the inverse

            IndexType bmpPos = 3 * ( i * width + j );
            IndexType imgPos = 3 * ( ( height - 1 - i ) * width + j );

            tmpData[ bmpPos   ] = static_cast<unsigned char>( rPixelData[imgPos + 2] );
            tmpData[ bmpPos + 1 ] = static_cast<unsigned char>( rPixelData[imgPos + 1] );
            tmpData[ bmpPos + 2 ] = static_cast<unsigned char>( rPixelData[imgPos] );
        }
    }

    for ( IndexType i = 0; i < height ; i++ )
    {
        fwrite( &tmpData[3 * width * i], 3, width, f );
    }

    fclose( f );
}

/* ------------------------------------------------------------------------------------ */
/*   Read PNG file                                                                      */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void ImageIO::readPNG( HArray<ValueType>& imageData, common::Grid3D& imageSize, const std::string& inputFileName )
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

    /* initialize stuff */
    png_structp png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
    SCAI_ASSERT_ERROR( png_ptr, "png_create_read_struct failed" )
    png_infop info_ptr = png_create_info_struct( png_ptr );
    SCAI_ASSERT_ERROR( info_ptr, "png_create_info_struct failed" )
    SCAI_ASSERT_ERROR( !setjmp( png_jmpbuf( png_ptr ) ), " FAIL" )
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
    SCAI_ASSERT_EQ_ERROR( 1, number_of_passes, "multiple passes not supported" )

    png_read_update_info( png_ptr, info_ptr );

    /* read file */

    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        COMMON_THROWEXCEPTION( "[read_png_file] Error during read_image" );
    }

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
/*   Read image file, make choice by suffix                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void ImageIO::writePNG( const HArray<ValueType>& data, const common::Grid3D& grid, const std::string& outputFileName )
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
    SCAI_ASSERT_ERROR( !setjmp( png_jmpbuf( png_ptr ) ), "failed" )
    png_init_io( png_ptr, fp );

    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        COMMON_THROWEXCEPTION( "[write png file error]" )
    }

    /* write header */
    png_byte bit_depth = 8;
    png_byte color_type = PNG_COLOR_TYPE_RGBA;
    png_set_IHDR( png_ptr, info_ptr, width, height,
                  bit_depth, color_type, PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );
    png_write_info( png_ptr, info_ptr );

    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        COMMON_THROWEXCEPTION( "[write png file error]" )
    }

    png_write_image( png_ptr, row_pointers );

    if ( setjmp( png_jmpbuf( png_ptr ) ) )
    {
        COMMON_THROWEXCEPTION( "[write png file error]" )
    }

    png_write_end( png_ptr, NULL );
    fclose( fp );
}

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
        readPNG( data, grid, inputFileName );
        SCAI_LOG_INFO( logger, "read PNG file " << inputFileName << ", data = " << data << ", grid = " << grid )
        imageData.swap( data, grid );
    }
    else if ( suffix == ".bmp" )
    {
        LArray<ValueType> data;
        readBMP( data, grid, inputFileName );
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
        writePNG( pixelData, pixelGrid, outputFileName );
    }
    else if ( suffix == ".bmp" )
    {
        writeBMP( pixelData, pixelGrid, outputFileName );
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
