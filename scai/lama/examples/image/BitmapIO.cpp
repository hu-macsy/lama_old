/**
 * @file BitmapIO.cpp
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

#include <scai/lama/examples/image/BitmapIO.hpp>

#include<png.h>
#include<fstream>

namespace scai
{

using namespace utilskernel;
using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( BitmapIO::logger, "ImageIO.Bitmap" )

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
/*   Read Bitmap File                                                                   */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void BitmapIO::readImpl( HArray<ValueType>& data, common::Grid3D& grid, const std::string& inputFileName )
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
    SCAI_LOG_INFO( logger, "infoHeader, biWidth = " << infoHeader.biWidth )
    SCAI_LOG_INFO( logger, "infoHeader, biHeight = " << infoHeader.biHeight )
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

    if ( infoHeader.biCompression == 3 )
    {
        // file.read( colorMasks, sizeof( colorMasks ) );
    }

    if ( infoHeader.biBitCount == 24 )
    {
        grid = common::Grid3D( height, width, 3 );

        // due to alignment: imageSize might be slightly higher

        SCAI_ASSERT_GE_ERROR( imageSize, grid.size(), "size mismatch" )
        common::scoped_array<unsigned char> tmp( new unsigned char[ imageSize ] );
        IndexType nBytes = fread( tmp.get(), 1, imageSize, file );
        SCAI_ASSERT_EQ_ERROR( imageSize, nBytes, "could not read all expected image pixels" )

        // Image lines are always stored as a multiple of 4

        IndexType widthAligned = ( ( width + 3 ) / 4 ) * 4;
        SCAI_ASSERT_EQ_ERROR( imageSize, widthAligned * height * 3, "unexpected image size, widthAligned = " << widthAligned )

        {
            WriteOnlyAccess<ValueType> wData( data, grid.size() );

            for ( IndexType i = 0; i < height; ++i )
            {
                for ( IndexType j = 0; j < width; ++j )
                {
                    // BMP has data from top to bottom, but we need the inverse
                    IndexType bmpPos = 3 * ( i * ( widthAligned )  + j );
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
    else if ( infoHeader.biBitCount == 8 )
    {
        grid = common::Grid3D( height, width, 3 );

        // indexed colors, read in color table

        unsigned char ColorMasks [256][4];

        IndexType nBytes = fread( ColorMasks, 1, sizeof( ColorMasks) , file );
        SCAI_LOG_ERROR( logger, "read color mask, nbytes = " << nBytes )
        common::scoped_array<unsigned char> tmp( new unsigned char[ width * height ] );
        nBytes = fread( tmp.get(), 1, width * height, file );
        SCAI_LOG_ERROR( logger, "read bitmap data, nbytes = " << nBytes )
        SCAI_LOG_ERROR( logger, "allocate data, size = " << grid.size() )

        WriteOnlyAccess<ValueType> wData( data, grid.size() );
        {
            for ( IndexType i = 0; i < height; ++i )
            {
                for ( IndexType j = 0; j < width; ++j )
                {
                    // BMP has data from top to bottom, but we need the inverse
                    IndexType bmpPos = i * width + j;
                    unsigned char colorIndex = tmp[ bmpPos ];
                    IndexType imgPos = 3 * ( ( height - 1 - i ) * width + j );
                    wData[ imgPos   ] = static_cast<ValueType>( ColorMasks[colorIndex][2]);
                    wData[ imgPos + 1 ] = static_cast<ValueType>( ColorMasks[colorIndex][1] );
                    wData[ imgPos + 2 ] = static_cast<ValueType>( ColorMasks[colorIndex][0] );
                }
            }
        }
        SCAI_LOG_ERROR( logger, "set full colo values" )
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Unsupported biBitCount = " << infoHeader.biBitCount )
    }
}

template<typename ValueType>
void BitmapIO::writeImpl( const HArray<ValueType>& data, const common::Grid3D& grid, const std::string& outputFileName )
{
    FILE* f = fopen( outputFileName.c_str(), "wb" );

    IndexType height  = grid.size( 0 );
    IndexType width   = grid.size( 1 );
    IndexType nColor  = grid.size( 2 );

    IndexType alignSize    = 4;
    IndexType widthAligned = ( ( width + alignSize ) / alignSize ) * alignSize;

    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;

    IndexType headerBytes = sizeof( fileHeader ) + sizeof( infoHeader );
    IndexType imageBytes  = widthAligned * height * nColor;

    fileHeader.bfType     = 0x4D42;
    fileHeader.bfSize     = imageBytes + headerBytes;
    fileHeader.bfReserved = 0;
    fileHeader.bfOffBits  = headerBytes;

    infoHeader.biSize = sizeof( infoHeader );
    infoHeader.biWidth = width;
    infoHeader.biHeight = height;
    infoHeader.biPlanes = 1;
    infoHeader.biBitCount = 24;
    infoHeader.biCompression = 0;
    infoHeader.biSizeImage = imageBytes;
    infoHeader.biXPelsPerMeter = 13780;
    infoHeader.biYPelsPerMeter = 13780;
    infoHeader.biClrUsed = 0;
    infoHeader.biClrImportant = 0;

    fwrite( &fileHeader, 1, sizeof( fileHeader ), f );
    fwrite( &infoHeader, 1, sizeof( infoHeader ), f );

    common::scoped_array<unsigned char> tmpData( new unsigned char[ imageBytes ] );

    ReadAccess<ValueType> rPixelData( data );

    for ( IndexType i = 0; i < height; ++i )
    {
        for ( IndexType j = 0; j < width; ++j )
        {
            // BMP has data from top to bottom, but we need the inverse

            IndexType bmpPos = 3 * ( i * widthAligned + j );
            IndexType imgPos = 3 * ( ( height - 1 - i ) * width + j );

            tmpData[ bmpPos   ] = static_cast<unsigned char>( rPixelData[imgPos + 2] );
            tmpData[ bmpPos + 1 ] = static_cast<unsigned char>( rPixelData[imgPos + 1] );
            tmpData[ bmpPos + 2 ] = static_cast<unsigned char>( rPixelData[imgPos] );
        }
    }

    IndexType nBytes = fwrite( tmpData.get(), 1, imageBytes, f );

    SCAI_ASSERT_EQ_ERROR( nBytes, imageBytes, "could not write all image data" )

    fclose( f );
}

/* ------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------ */

void BitmapIO::read( _HArray& data, common::Grid3D& grid, const std::string& inputFileName )
{
    if ( data.getValueType() == common::scalar::FLOAT )
    {
        HArray<float>& dataF = reinterpret_cast<HArray<float>&>( data );
        readImpl( dataF, grid, inputFileName );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported data type " << data.getValueType() );
    }
}

void BitmapIO::write( const _HArray& data, const common::Grid3D& grid, const std::string& inputFileName )
{
    if ( data.getValueType() == common::scalar::FLOAT )
    {
        const HArray<float>& dataF = reinterpret_cast<const HArray<float>&>( data );
        writeImpl( dataF, grid, inputFileName );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported data type " << data.getValueType() );
    }
}

// instantiate methods

// template
// void BitmapIO::readImpl( HArray<float>&, common::Grid3D& grid, const std::string& );


} /* end namespace lama */

} /* end namespace scai */
