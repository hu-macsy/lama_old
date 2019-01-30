/**
 * @file BitmapIO.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Implementation of routines to read/write image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/GridVector.hpp>

#include <scai/lama/io/BitmapIO.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/common/exception/UnsupportedException.hpp>

#include <fstream>

#define MAT_SUFFIX ".bmp" 

namespace scai
{

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( BitmapIO::logger, "ImageIO.Bitmap" )

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* BitmapIO::create()
{
    return new BitmapIO();
}

std::string BitmapIO::createValue()
{
    return MAT_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool BitmapIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == FileMode::FORMATTED )
    {
        return false;
    }

    return true;
}

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

/* --------------------------------------------------------------------------------- */

void BitmapIO::open( const char* fileName, const char* fileMode )
{
    if ( strcmp( fileMode, "w" ) == 0 )
    {
        mFile = fopen( fileName, "wb" );
    }
    else if ( strcmp( fileMode, "r" ) == 0 )
    {
        mFile = fopen( fileName, "rb" );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported file mode for Bitmap file: " << fileMode )
    }

    SCAI_ASSERT( mFile, "Could not open bitmap file " << fileName << ", mode = " << fileMode )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::close()
{
    fclose( mFile );
}

/* ------------------------------------------------------------------------------------ */
/*   Read Bitmap File                                                                   */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void BitmapIO::readGridImpl( HArray<ValueType>& data, common::Grid& grid )
{
    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;
    // COLORMASKS colorMasks;
    // COLORTABLE colorTable;
    size_t nRecords = fread( &fileHeader, sizeof( fileHeader ), 1, mFile );

    SCAI_ASSERT_EQ_ERROR( 1, nRecords, "failed to read header" )
    SCAI_ASSERT_EQ_ERROR( fileHeader.bfType, 0x4D42, "not bitmap file" )


    SCAI_LOG_INFO( logger, "fileHeader, bfSize = " << fileHeader.bfSize
                           << ", bfReserved = " << fileHeader.bfReserved 
                           << ", bfOffBits = " << fileHeader.bfOffBits )


    nRecords = fread( &infoHeader, sizeof( infoHeader ), 1, mFile );

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

    if ( imageSize == 0 )
    {
        imageSize = fileHeader.bfSize - fileHeader.bfOffBits;
        SCAI_LOG_INFO( logger, "Take image size = " << imageSize << " from fileHeader" )
    }

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
        std::unique_ptr<unsigned char[]> tmp( new unsigned char[ imageSize ] );
        IndexType nBytes = fread( tmp.get(), 1, imageSize, mFile );
        SCAI_ASSERT_EQ_ERROR( imageSize, nBytes, "could not read all expected image pixels" )

        // Image lines are always stored as a multiple of 4

        IndexType widthAligned = ( ( width + 3 ) / 4 ) * 4;

        // SCAI_ASSERT_EQ_ERROR( imageSize, widthAligned * height * 3, "unexpected image size, widthAligned = " << widthAligned )

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

        IndexType nBytes = fread( ColorMasks, 1, sizeof( ColorMasks), mFile );
        SCAI_LOG_INFO( logger, "read color mask, nbytes = " << nBytes )
        std::unique_ptr<unsigned char[]> tmp( new unsigned char[ width * height ] );
        nBytes = fread( tmp.get(), 1, width * height, mFile );
        SCAI_LOG_INFO( logger, "read bitmap data, nbytes = " << nBytes )
        SCAI_LOG_INFO( logger, "allocate data, size = " << grid.size() )

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
        SCAI_LOG_INFO( logger, "set full colo values" )
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Unsupported biBitCount = " << infoHeader.biBitCount )
    }
}

template<typename ValueType>
void BitmapIO::writeGridImpl( const HArray<ValueType>& data, const common::Grid& grid )
{
    SCAI_ASSERT_EQ_ERROR( 3, grid.nDims(), "image data must be 3-dimensional grid data" )

    IndexType height  = grid.size( 0 );
    IndexType width   = grid.size( 1 );
    IndexType nColor  = grid.size( 2 );

    SCAI_ASSERT_EQ_ERROR( 3, nColor, "only RGB images are supported" )

    IndexType alignSize    = 4;
    IndexType widthAligned = ( ( width * nColor  + alignSize - 1 ) / alignSize ) * alignSize;

    SCAI_LOG_INFO( logger, "row size = " << width * nColor << ", aligned = " << widthAligned )

    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;

    IndexType headerBytes = sizeof( fileHeader ) + sizeof( infoHeader );
    IndexType imageBytes  = widthAligned * height;

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

    fwrite( &fileHeader, 1, sizeof( fileHeader ), mFile );
    fwrite( &infoHeader, 1, sizeof( infoHeader ), mFile );

    std::unique_ptr<unsigned char[]> tmpData( new unsigned char[ imageBytes ] );

    ReadAccess<ValueType> rPixelData( data );

    for ( IndexType i = 0; i < height; ++i )
    {
        for ( IndexType j = 0; j < width; ++j )
        {
            // BMP has data from top to bottom, but we need the inverse

            IndexType bmpPos = i * widthAligned + nColor * j ;
            IndexType imgPos = nColor * ( ( height - 1 - i ) * width + j );

            tmpData[ bmpPos   ] = static_cast<unsigned char>( rPixelData[imgPos + 2] );
            tmpData[ bmpPos + 1 ] = static_cast<unsigned char>( rPixelData[imgPos + 1] );
            tmpData[ bmpPos + 2 ] = static_cast<unsigned char>( rPixelData[imgPos] );
        }
    }

    IndexType nBytes = fwrite( tmpData.get(), 1, imageBytes, mFile );

    SCAI_ASSERT_EQ_ERROR( nBytes, imageBytes, "could not write all image data" )
}

/* ------------------------------------------------------------------------------------ */

void BitmapIO::readGridArray( _HArray& data, common::Grid& grid )
{
    IOWrapper<BitmapIO, SCAI_TYPELIST( float, double )>::readGridImpl( *this, data, grid );
}

void BitmapIO::writeGridArray( const _HArray& data, const common::Grid& grid )
{
    IOWrapper<BitmapIO, SCAI_TYPELIST( float, double )>::writeGridImpl( *this, data, grid );
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::writeStorage( const _MatrixStorage& storage )
{
    // todo: write dense storage as bitmap, maybe greyscale 

    SCAI_THROWEXCEPTION( common::UnsupportedException, "write storage " << storage )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::readStorage( _MatrixStorage& storage )
{
    storage.clear();

    SCAI_THROWEXCEPTION( common::UnsupportedException, "Unsupported for bitmap file: read storage" )
}

/* --------------------------------------------------------------------------------- */

std::string BitmapIO::getMatrixFileSuffix() const
{
    return createValue();
}

/* --------------------------------------------------------------------------------- */

std::string BitmapIO::getVectorFileSuffix() const
{
    return createValue();
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::writeArray( const hmemo::_HArray& array )
{
    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: write array " << array )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::writeSparse( 
    const IndexType n, 
    const hmemo::HArray<IndexType>& indexes, 
    const hmemo::_HArray& values )
{
    SCAI_THROWEXCEPTION( common::UnsupportedException, 
                         "Unsupported for bitmap file: write spare array ( n = " << n 
                         << ", indexes = " << indexes << ", values = " << values << " )" )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    numRows    = 0;
    numColumns = 0;
    numValues  = 0;

    COMMON_THROWEXCEPTION( "Unsupported for bitmap file: read storage info" )
}

void BitmapIO::getArrayInfo( IndexType& size )
{
    size = 0;

    SCAI_THROWEXCEPTION( common::UnsupportedException, "Unsupported for bitmap file: read array info" )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::readArray( hmemo::_HArray& array )
{
    array.clear();

    SCAI_THROWEXCEPTION( common::UnsupportedException, "Unsupported for bitmap file: read array" )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    size = 0;
    indexes.clear();
    values.clear();

    SCAI_THROWEXCEPTION( common::UnsupportedException, "Unsupported for bitmap file: read sparse array" )
}

/* --------------------------------------------------------------------------------- */

void BitmapIO::writeAt( std::ostream& stream ) const
{
    stream << "BitmapIO ( suffix = " << MAT_SUFFIX << ", ";
    stream << ", only image data )";
}

/* --------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
