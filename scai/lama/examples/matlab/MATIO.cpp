/**
 * @file MATIO.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of methods for FileIO class MATIO
 * @author Thomas Brandes
 * @date 10.06.2016
 */


#include "MATIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/io/IOStream.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/exception/IOException.hpp>

#include <sstream>
#include <zlib.h>

#define MAT_SUFFIX ".mat"

using namespace std;

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* MATIO::create()
{
    return new MATIO();
}

std::string MATIO::createValue()
{
    return MAT_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool MATIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void MATIO::writeAt( std::ostream& stream ) const
{
    stream << "MATIO ( suffix = " << MAT_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MATIO::logger, "FileIO.MATIO" )

/* --------------------------------------------------------------------------------- */

static void readMATFileHeader(
    int& version,
    IOStream::Endian& indicator,
    IOStream& inFile )
{
    // MATFile header has 128 bytes

    char header[128];

    inFile.read( header, 128 );

    header[116] = '\0';

    const short int* versionPtr   = reinterpret_cast<short int*>( &header[124] );
    // const short int* indicatorPtr = reinterpret_cast<short int*>( &header[126] );

    std::cout << "Header txt = " << header << std::endl;

    indicator = IOStream::MACHINE_ENDIAN;
    version = *versionPtr;

    std::cout << "Version = " << version << std::endl;
}

/* --------------------------------------------------------------------------------- */

static const char* readDataElementHeader( uint32_t& dataType, uint32_t& nBytes, uint32_t& wBytes, const char buffer[] )
{
    const uint32_t* data = reinterpret_cast<const uint32_t*>( buffer );

    bool small = false;

    dataType = data[0];

    if ( ( dataType >> 16 ) != 0 )
    {
        small    = true;
        nBytes   = dataType >> 16;
        dataType = ( dataType << 16 ) >> 16;
        std::cout << "small element: nBytes = " << nBytes << ", dataType = " << dataType << std::endl;
    }
    else
    {
        nBytes   = data[1];
    }

    if ( small )
    {
        wBytes = 8;
        return buffer + 4;   // the next 4 bytes are just the data
    }
    else
    {
        wBytes = 8 + nBytes + ( 8 - nBytes ) % 8;
        return buffer + 8;
    }
}

/* --------------------------------------------------------------------------------- */

static void uncompress( char out[], int out_len, const char in[], int len, int flush )
{
    z_stream infstream;

    infstream.zalloc = Z_NULL;
    infstream.zfree  = Z_NULL;
    infstream.opaque = Z_NULL;

    int ok = inflateInit( &infstream );

    SCAI_ASSERT_EQ_ERROR( ok, Z_OK, "Unable to inflateInit" )

    // inflate the variable head

    infstream.avail_in = len;
    infstream.next_in = ( unsigned char* ) in;
    infstream.avail_out = out_len;
    infstream.next_out = ( unsigned char* ) out;

    ok = inflate( &infstream, flush );

    if ( flush == Z_NO_FLUSH )
    {
        SCAI_ASSERT_EQ_ERROR( ok, Z_OK, "Unable to inflate" )
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( ok, Z_STREAM_END, "Unable to inflate" )
    }

    // get the headers

    // inflate(&infstream, Z_FINISH);
    inflateEnd( &infstream );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MATIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    IOStream outFile( fileName, std::ios::out );

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( array, precData );
}

/* --------------------------------------------------------------------------------- */

void MATIO::readArrayInfo( IndexType&, const std::string& )
{
    COMMON_THROWEXCEPTION( "not available" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MATIO::readArrayImpl(
    hmemo::HArray<ValueType>& ,
    const std::string& ,
    const IndexType ,
    const IndexType )
{
    COMMON_THROWEXCEPTION( "not available" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MATIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    COOStorage<ValueType> coo( storage );

    HArray<IndexType> cooIA = coo.getIA();
    HArray<IndexType> cooJA = coo.getJA();
    HArray<ValueType> cooValues = coo.getValues();

    IOStream outFile( fileName, std::ios::out );

    int precIndex = 0;
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    if ( mScalarTypeData == common::scalar::PATTERN )
    {
        outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex );
    }
    else
    {
        outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex, cooValues, precData );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MATIO::readData(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>* values,
    const IndexType nnz,
    const std::string& fileName )
{
    LArray<double> dIA;
    LArray<double> dJA;

    IOStream inFile( fileName, std::ios::in );

    if ( values != NULL )
    {
        inFile.readFormatted( dIA, dJA, *values, nnz );
    }
    else
    {
        inFile.readFormatted( dIA, dJA, nnz );
    }

    ContextPtr ctx = Context::getHostPtr();

    HArrayUtils::setArrayImpl( ia, dIA );  // conversion from double to IndexType
    HArrayUtils::setArrayImpl( ja, dJA );  // conversion from double to IndexType

    IndexType minRowIndex = HArrayUtils::reduce( ia, binary::MIN );

    if ( minRowIndex == 0 )
    {
        // okay, seems that indexing start with 0
    }
    else if ( minRowIndex == 1 )
    {
        // offset base = 1, convert it to 0

        HArrayUtils::setScalar( ia, IndexType( 1 ), binary::SUB );
        HArrayUtils::setScalar( ja, IndexType( 1 ), binary::SUB );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ERROR reading file " << fileName << ": minimal row index " << minRowIndex << " is illegal" )
    }
}

/* --------------------------------------------------------------------------------- */

void MATIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName )
{
    numRows    = 0;
    numColumns = 0;
    numValues  = 0;

    COMMON_THROWEXCEPTION( "not read " << fileName )
}

/* --------------------------------------------------------------------------------- */

static void getStorage( _MatrixStorage& storage, const char* dataElementPtr, uint32_t dataType, uint32_t nBytes )
{
    std::cout << "getStorage, dataType = " << dataType << ", nBytes = " << nBytes << std::endl;

    if ( dataType == 14 )
    {
        // read subelements of miMatrix

        uint32_t offset = 8;

        uint32_t subType;
        uint32_t subNBytes;
        uint32_t subWBytes;

        readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );

        SCAI_ASSERT_EQ_ERROR( subNBytes, 8, "ArrayFlags must be 8 bytes" )
        SCAI_ASSERT_EQ_ERROR( subType, 6, "ArrayFlags type must be miUINIT32" )

        uint8_t flags = dataElementPtr[ offset + 9 ];
        uint8_t cbyte = dataElementPtr[ offset + 8 ];

        bool isComplex = flags & ( 1 << 3 );

        const int* nzPtr = reinterpret_cast<const int*>( dataElementPtr + offset + 12 );

        uint32_t nzmax = *nzPtr;

        std::cout << "array, flags = " << ( int ) flags << ", class = " << ( int ) cbyte << ", nzmax = " << nzmax << std::endl;

        offset += 16;

        readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );

        std::cout << "array, dims, type = " << subType << ", #bytes = " << subNBytes << ", #offset = " << subWBytes << std::endl;

        const int* dimPtr = reinterpret_cast<const int*>( dataElementPtr + offset + 8 );

        std::cout << "array, dims = " << dimPtr[0] << " x " << dimPtr[1] << std::endl;

        offset += 16;

        // now read the name

        const char* name = readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );

        std::cout << "array, name, type = " << subType << ", #bytes = " << subNBytes << std::endl;

        common::scoped_array<char> nameData( new char[ subNBytes + 1 ] );

        for ( uint32_t i = 0; i < subNBytes; ++i )
        {
            nameData[i] = name[i];
        }

        nameData[subNBytes] = '\0';

        std::cout << "array name = " << nameData.get() << std::endl;

        offset += subWBytes;

        if ( cbyte == 5 )
        {
            LArray<IndexType> ia;
            LArray<IndexType> ja;
            LArray<double> values;

            const void* iaPtr = readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
            std::cout << "sparse row indexes, type = " << subType << ", #bytes = " << subNBytes << std::endl;
            SCAI_ASSERT_EQ_ERROR( subType, 5, "row index type mismatch" )
            SCAI_ASSERT_EQ_ERROR( subNBytes, nzmax * 4, "row index size mismatch" )
            {
                hmemo::WriteOnlyAccess<IndexType> wIA( ia, nzmax );
                ::memcpy( wIA.get(), iaPtr, subNBytes );
            }
            offset += subWBytes;
            const void* jaPtr = readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
            std::cout << "sparse col indexes, type = " << subType << ", #bytes = " << subNBytes << std::endl;
            SCAI_ASSERT_EQ_ERROR( subType, 5, "col index type mismatch" )
            SCAI_ASSERT_EQ_ERROR( subNBytes, static_cast<uint32_t>( dimPtr[1] + 1 ) * 4, "col index size mismatch" )
            {
                // ja is an offset array that is translated to nzmax column indexes
                // e.g. 0, 2, 5, 7, 8  -> 0 0 1 1 1 2 2 3
                const IndexType* offsetJA = reinterpret_cast<const IndexType*>( jaPtr );
                static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::offsets2ia> offsets2ia;
                hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
                offsets2ia.getSupportedContext( loc );
                hmemo::WriteOnlyAccess<IndexType> wJA( ja, loc, nzmax );
                SCAI_CONTEXT_ACCESS( loc )
                offsets2ia[loc]( wJA.get(), nzmax, offsetJA, dimPtr[1], 0 );
            }

            offset += subWBytes;
            const void* valuesPtr = readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
            std::cout << "sparse real values = " << subType << ", #bytes = " << subNBytes << std::endl;
            SCAI_ASSERT_EQ_ERROR( subType, 9, "sparse values type mismatch" )
            SCAI_ASSERT_EQ_ERROR( subNBytes, nzmax * 8, "value size mismatch" )
            {
                hmemo::WriteOnlyAccess<double> wValues( values, nzmax );
                ::memcpy( wValues.get(), valuesPtr, subNBytes );
            }
            offset += subWBytes;

            if ( isComplex )
            {
                readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
                std::cout << "sparse imag values = " << subType << ", #bytes = " << subNBytes << std::endl;
                SCAI_ASSERT_EQ_ERROR( subType, 9, "sparse imag values type mismatch" )
                SCAI_ASSERT_EQ_ERROR( subNBytes, nzmax * 8, "sparse imag value size mismatch" )
                offset += subWBytes;
            }

            COOStorage<double> cooStorage( dimPtr[0], dimPtr[1] );
            cooStorage.swap( ia, ja, values );
            storage = cooStorage;
        }
        else
        {
            readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
            std::cout << "array values, type = " << subType << ", #bytes = " << subNBytes << std::endl;
            offset += subWBytes;

            if ( isComplex )
            {
                readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );
                std::cout << "array complex values, type = " << subType << ", #bytes = " << subNBytes << std::endl;
                offset += subWBytes;
            }
        }

        std::cout << "Final offset = " << offset << ", len variable = " << nBytes << endl;

        SCAI_ASSERT_EQ_ERROR( offset, nBytes, "mismatch read bytes and size bytes, maybe COMPLEX" )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MATIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& matrixFileName,
    const IndexType ,
    const IndexType )
{
    IOStream inFile( matrixFileName, std::ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;
    readMATFileHeader( version, endian, inFile );

    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    common::scoped_array<char> dataElement;

    const char* dataElementPtr;

    while ( true )
    {
        char buffer[8];

        inFile.read( buffer, 8 );

        if ( inFile.fail() )
        {
            std::cout << "No more header" << std::endl;
            break;
        }

        readDataElementHeader( dataType, nBytes, wBytes, buffer );

        std::cout << std::endl;
        std::cout << "Next data element in MAT file" << std::endl;
        std::cout << "=============================" << std::endl;
        std::cout << "dataType = " << dataType << std::endl;
        std::cout << "nBytes = " << nBytes << std::endl;

        if ( dataType == 15 )
        {
            common::scoped_array<char> compressedData( new char[nBytes] );

            inFile.read( compressedData.get(), nBytes );

            uncompress( buffer, 8, compressedData.get(), nBytes, Z_NO_FLUSH );

            uint32_t uncompressedNBytes;
            uint32_t uncompressedWBytes;

            readDataElementHeader( dataType, uncompressedNBytes, uncompressedWBytes, buffer );

            std::cout << "uncompressed dataType = " << dataType << std::endl;
            std::cout << "uncompressed nBytes = " << uncompressedNBytes << std::endl;

            dataElement.reset( new char[uncompressedNBytes + 8] );

            uncompress( dataElement.get(), uncompressedNBytes + 8, compressedData.get(), nBytes, Z_FINISH );

            dataElementPtr = dataElement.get();

            nBytes = uncompressedNBytes + 8;
        }
        else
        {
            dataElement.reset( new char[ nBytes] );
            inFile.read( dataElement.get(), nBytes );
            dataElementPtr = dataElement.get();
        }

        getStorage( storage, dataElementPtr, dataType, nBytes );

        if ( inFile.eof() )
        {
            break;
        }
        else
        {
            std::cout << "there is another data element." << std::endl;
        }
    }

    std::cout << "storage = " << storage << std::endl;
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
