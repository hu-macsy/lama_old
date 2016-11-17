/**
 * @file MatlabIO.cpp
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
 * @brief Implementation of methods for FileIO class MatlabIO
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include "MatlabIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
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

enum {
    MAT_INT8       = 1,
    MAT_UINT8      = 2,
    MAT_INT16      = 3,
    MAT_UINT16     = 4,
    MAT_INT32      = 5,
    MAT_UINT32     = 6,
    MAT_FLOAT      = 7,
    MAT_DOUBLE     = 9,
    MAT_INT64      = 12,
    MAT_UINT64     = 13,
    MAT_MATRIX     = 14,
    MAT_COMPRESSED = 15,
    MAT_UTF8       = 16,
    MAT_UTF16      = 17,
    MAT_UTF32      = 18
};

enum {
    MAT_CELL_CLASS     = 1,
    MAT_STRUCT_CLASS   = 2,
    MAT_OBJECT_CLASS   = 3,
    MAT_CHAR_CLASS     = 4,
    MAT_SPARSE_CLASS   = 5,
    MAT_DOUBLE_CLASS   = 6,
    MAT_FLOAT_CLASS    = 7,
    MAT_INT8_CLASS     = 8,
    MAT_UINT8_CLASS    = 9,
    MAT_INT16_CLASS    = 10,
    MAT_UINT16_CLASS   = 11,
    MAT_INT32_CLASS    = 12,
    MAT_UINT32_CLASS   = 13,
    MAT_INT64_CLASS    = 14,
    MAT_UINT64_CLASS   = 15
};

/* --------------------------------------------------------------------------------- */

static common::scalar::ScalarType matlabType2ScalarType( uint32_t dataType )
{
    switch( dataType )
    {
        case  MAT_INT8   : return common::scalar::CHAR;          // INT8 same as CHAR
        case  MAT_UINT8  : return common::scalar::UNKNOWN;       // UINT8 not supported yet
        case  MAT_INT16  : return common::scalar::UNKNOWN;       // INT16 not supported yet
        case  MAT_UINT16 : return common::scalar::UNKNOWN;       // UINT16 not supported yet
        case  MAT_INT32  : return common::scalar::INT;           // INT32
        case  MAT_UINT32 : return common::scalar::UNSIGNED_INT;  // UINT32 not supported yet
        case  MAT_FLOAT  : return common::scalar::FLOAT;         // single precision
        case  MAT_DOUBLE : return common::scalar::DOUBLE;        // double precision
        case  MAT_INT64  : return common::scalar::LONG;          // INT64
        case  MAT_UINT64 : return common::scalar::UNSIGNED_LONG; // UINT64
        default          : return common::scalar::UNKNOWN;       // 
    }
}

/* --------------------------------------------------------------------------------- */

static common::scalar::ScalarType class2ScalarType( uint32_t dataType )
{
    switch( dataType )
    {
        case  MAT_INT8_CLASS   : return common::scalar::CHAR;          // INT8 same as CHAR
        case  MAT_UINT8_CLASS  : return common::scalar::UNKNOWN;       // UINT8 not supported yet
        case  MAT_INT16_CLASS  : return common::scalar::UNKNOWN;       // INT16 not supported yet
        case  MAT_UINT16_CLASS : return common::scalar::UNKNOWN;       // UINT16 not supported yet
        case  MAT_INT32_CLASS  : return common::scalar::INT;           // INT32
        case  MAT_UINT32_CLASS : return common::scalar::UNSIGNED_INT;  // UINT32 not supported yet
        case  MAT_FLOAT_CLASS  : return common::scalar::FLOAT;         // single precision
        case  MAT_DOUBLE_CLASS : return common::scalar::DOUBLE;        // double precision
        case  MAT_INT64_CLASS  : return common::scalar::LONG;          // INT64
        case  MAT_UINT64_CLASS : return common::scalar::UNSIGNED_LONG; // UINT64
        default                : return common::scalar::UNKNOWN;       // 
    }
}

static uint8_t scalarType2Class( common::scalar::ScalarType stype ) 
{
    switch( stype )
    {
        case  common::scalar::CHAR           : return MAT_INT8_CLASS;
        case  common::scalar::INT            : return MAT_INT32_CLASS;
        case  common::scalar::UNSIGNED_INT   : return MAT_UINT32_CLASS;
        case  common::scalar::LONG           : return MAT_INT64_CLASS;
        case  common::scalar::UNSIGNED_LONG  : return MAT_UINT64_CLASS;
        case  common::scalar::FLOAT          : return MAT_FLOAT_CLASS;
        case  common::scalar::DOUBLE         : return MAT_DOUBLE_CLASS;
        case  common::scalar::COMPLEX        : return MAT_FLOAT_CLASS;
        case  common::scalar::DOUBLE_COMPLEX : return MAT_DOUBLE_CLASS;

        default                            : COMMON_THROWEXCEPTION( "no Matlab class for " << stype )
    }

    return 0;
}

template<typename ValueType> 
uint32_t matlabType()
{
    COMMON_THROWEXCEPTION( "unsupported type for Matlab: " << typeid( ValueType ).name() )
    return 0;
}

template<>
uint32_t matlabType<char>()
{
    return MAT_INT8;
}

template<>
uint32_t matlabType<int>()
{
    return MAT_INT32;
}

template<>
uint32_t matlabType<unsigned int>()
{
    return MAT_UINT32;
}

template<>
uint32_t matlabType<long>()
{
    return MAT_INT64;
}

template<>
uint32_t matlabType<unsigned long>()
{
    return MAT_UINT64;
}

template<>
uint32_t matlabType<float>()
{
    return MAT_FLOAT;
}

template<>
uint32_t matlabType<double>()
{
    return MAT_DOUBLE;
}

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* MatlabIO::create()
{
    return new MatlabIO();
}

std::string MatlabIO::createValue()
{
    return MAT_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool MatlabIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeAt( std::ostream& stream ) const
{
    stream << "MatlabIO ( suffix = " << MAT_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatlabIO::logger, "FileIO.MatlabIO" )

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
        std::cout << "data element: nBytes = " << nBytes << " / " << wBytes << ", dataType = " << dataType << std::endl;
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

static void readDataElement( IOStream& inFile, common::scoped_array<char>& dataElement )
{
    char buffer[8];

    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    inFile.read( buffer, 8 );

    if ( inFile.fail() )
    {
        COMMON_THROWEXCEPTION( "No element found in input file" );
    }

    readDataElementHeader( dataType, nBytes, wBytes, buffer );

    if ( dataType != MAT_COMPRESSED )
    {
        // uncompressed data, copy header and the other data

        std::cout << "uncompressed elem, width = " << wBytes << ", #bytes = " << nBytes << endl;

        dataElement.reset( new char[ wBytes ] );
        memcpy( dataElement.get(), buffer, 8 );
        inFile.read( dataElement.get() + 8, wBytes - 8 );
    }
    else
    {
        common::scoped_array<char> compressedData( new char[nBytes] );

        inFile.read( compressedData.get(), nBytes );

        uncompress( buffer, 8, compressedData.get(), nBytes, Z_NO_FLUSH );

        uint32_t uncompressedNBytes;
        uint32_t uncompressedWBytes;

        readDataElementHeader( dataType, uncompressedNBytes, uncompressedWBytes, buffer );

        std::cout << "uncompressed dataType = " << dataType << std::endl;
        std::cout << "uncompressed nBytes = " << uncompressedNBytes << std::endl;

        dataElement.reset( new char[uncompressedWBytes] );

        uncompress( dataElement.get(), uncompressedWBytes, compressedData.get(), nBytes, Z_FINISH );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType> 
static uint32_t getData( ValueType* data, uint32_t size, const char* buffer  )
{
    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    std::cout << "getData( " << common::TypeTraits<ValueType>::id() << ", size = " << size << std::endl;

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, buffer );

    SCAI_ASSERT_EQ_ERROR( nBytes, size * sizeof( ValueType), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( dataType, matlabType<ValueType>(), "type mismatch" )

    ::memcpy( data, dataPtr, nBytes );

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

static uint32_t getString( char* name, uint32_t nameSize, const char* buffer  )
{
    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, buffer );

    SCAI_ASSERT_EQ_ERROR( dataType, matlabType<char>(), "type mismatch" )
    SCAI_ASSERT_LT_ERROR( nBytes, nameSize - 1, "too long string" )

    ::memcpy( name, dataPtr, nBytes );

    name[nBytes] = '\0';   // finalize string

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

template<typename ArrayType, typename DataType>
static void readMATArrayImpl( hmemo::HArray<ArrayType>& array, const char* data, IndexType nBytes )
{
    IndexType elemSize  = sizeof( DataType );
    IndexType arraySize = nBytes / elemSize;

    SCAI_ASSERT_EQ_ERROR( elemSize * arraySize, nBytes, "Size mismatch, elemSize = " << elemSize << ", arraySize = " << arraySize )

    if ( typeid( ArrayType ) == typeid( DataType ) )
    {
        std::cout << "readMatArrayImpl, in place, type = " << common::TypeTraits<ArrayType>::id() 
                  << ", arraySize = " << arraySize << ", #bytes = " << nBytes << std::endl;

        // no temporary array required

        hmemo::WriteOnlyAccess<ArrayType> wData( array, arraySize );
        ::memcpy( wData.get(), data, nBytes );

        std::cout << "last = " << wData[arraySize-1] << std::endl;
    }
    else
    {
        std::cout << "readMatArrayImpl, convert, array type = " << common::TypeTraits<ArrayType>::id() 
                  << ", data type = " << typeid( DataType ).name() << ", array size = " << arraySize << std::endl;

        // temporary array and conversion required

        HArray<DataType> tmp;

        {
            hmemo::WriteOnlyAccess<DataType> wData( tmp, arraySize );
            ::memcpy( wData.get(), data, nBytes );
        }

        HArrayUtils::assign( array, tmp );  
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
static void readMATArray( hmemo::HArray<ValueType>& array, const char* data, const uint32_t mxClass, const uint32_t nbytes )
{
    switch ( mxClass )
    {
        case  MAT_DOUBLE_CLASS : readMATArrayImpl<ValueType, double>( array, data, nbytes ); break;
        case  MAT_FLOAT_CLASS  : readMATArrayImpl<ValueType, float>( array, data, nbytes ); break;
        case  MAT_INT8_CLASS   : readMATArrayImpl<ValueType, int8_t>( array, data, nbytes ); break;
        case  MAT_UINT8_CLASS  : readMATArrayImpl<ValueType, uint8_t>( array, data, nbytes ); break;
        case  MAT_INT16_CLASS  : readMATArrayImpl<ValueType, int16_t>( array, data, nbytes ); break;
        case  MAT_UINT16_CLASS : readMATArrayImpl<ValueType, uint16_t>( array, data, nbytes ); break;
        case  MAT_INT32_CLASS  : readMATArrayImpl<ValueType, int32_t>( array, data, nbytes ); break;
        case  MAT_UINT32_CLASS : readMATArrayImpl<ValueType, uint32_t>( array, data, nbytes ); break;
        case  MAT_INT64_CLASS  : readMATArrayImpl<ValueType, int64_t>( array, data, nbytes ); break;
        case  MAT_UINT64_CLASS : readMATArrayImpl<ValueType, uint64_t>( array, data, nbytes ); break;

        default : COMMON_THROWEXCEPTION( "mxClass = " << mxClass << " is unknown array class in Matlab file." )
    }
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readArrayInfo( IndexType& n, const std::string& arrayFileName )
{
    IOStream inFile( arrayFileName, std::ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    readMATFileHeader( version, endian, inFile );

    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    common::scoped_array<char> dataElement;

    readDataElement( inFile, dataElement );

    const char* elementPtr = readDataElementHeader( dataType, nBytes, wBytes, dataElement.get() );

    SCAI_ASSERT_EQ_ERROR( dataType, MAT_MATRIX, "can only read array as MATRIX - MATLAB array" )

    uint32_t header[2];
    int32_t  dims[2];

    uint32_t offset = getData( header, 2, elementPtr );
    offset += getData( dims, 2, elementPtr + offset );

    n = dims[0] * dims[1];
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& arrayFileName,
    const IndexType ,
    const IndexType )
{
    IOStream inFile( arrayFileName, std::ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    readMATFileHeader( version, endian, inFile );

    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    common::scoped_array<char> dataElement;

    readDataElement( inFile, dataElement );

    cout << "Read full data element from input file." << endl;

    const char* elementPtr = readDataElementHeader( dataType, nBytes, wBytes, dataElement.get() );

    SCAI_ASSERT_EQ_ERROR( dataType, MAT_MATRIX, "can only read array as MATRIX - MATLAB array" )

    uint32_t header[2];

    uint32_t offset = getData( header, 2, elementPtr );

    uint8_t* headerFlags = reinterpret_cast<uint8_t*>( header );

    // read flags and class of miMatrix

    uint8_t flags = headerFlags[ 1 ];
    uint8_t cbyte = headerFlags[ 0 ];

    // assert mxDOUBLE_CLASS

    bool isComplex = flags & ( 1 << 3 );

    uint32_t nzmax = header[1];

    std::cout << "array, flags = " << ( int ) flags << ", class = " << ( int ) cbyte 
              << ", nzmax = " << nzmax << ", isComplex = " << isComplex << std::endl;

    // SCAI_ASSERT_EQ_ERROR( static_cast<int>( cbyte ), 6, "not double precision array" )

    int dims[2];

    offset += getData( dims, 2, elementPtr + offset );

    std::cout << "array, dims = " << dims[0] << " x " << dims[1] << std::endl;

    // now read the name, maximal size is 128 characters

    char nameData[128];

    offset += getString( nameData, 128, elementPtr + offset );

    std::cout << "array name = " << nameData << std::endl;

    // now read the data

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, elementPtr + offset );
    offset += wBytes;

    std::cout << "array values, type = " << dataType << ", #bytes = " << nBytes << std::endl;

    SCAI_ASSERT_EQUAL( class2ScalarType( cbyte ), matlabType2ScalarType( dataType ), "type mismatch" )

    readMATArray( array, dataPtr, static_cast<uint32_t>( cbyte ), nBytes );

    SCAI_ASSERT_EQ_ERROR( array.size(), dims[0] * dims[1], "serious mismatch" )

    if ( isComplex )
    {
        dataPtr = readDataElementHeader( dataType, nBytes, wBytes, elementPtr + offset );
        offset += wBytes;

        std::cout << "array complex values, type = " << dataType << ", #bytes = " << nBytes << std::endl;

        if ( common::isComplex( array.getValueType() ) ) 
        {
            utilskernel::LArray<ValueType> imagValues;

            readMATArray( imagValues, dataPtr, static_cast<uint32_t>( cbyte ), nBytes );

            ValueType i = static_cast<ValueType>( ComplexDouble( 0, 1 ) );
            utilskernel::HArrayUtils::binaryOpScalar2( imagValues, imagValues, i, utilskernel::binary::MULT );
            utilskernel::HArrayUtils::binaryOp( array, array, imagValues, utilskernel::binary::ADD );
        }
        else
        {
            std::cout << "imaginary values of complex numbers ignored, array is not complex" << std::endl;
        }
    }
}

/* --------------------------------------------------------------------------------- */

static void writeMATFileHeader( IOStream& outFile )
{
    char header[128];

    memset( header, 0, 128 );

    sprintf( header, "This is the output of LAMA" );

    short int* versionPtr   = reinterpret_cast<short int*>( &header[124] );

    *versionPtr = 256;

    header[126] = 'I';
    header[127] = 'M';

    outFile.write( header, 128 );
}

/* --------------------------------------------------------------------------------- */

static void writeDataElementHeader( IOStream& outFile, const uint32_t dataType, const uint32_t nBytes )
{
    std::cout << "Write: DataElementHeader, type = " << dataType << " " << matlabType2ScalarType( dataType ) 
              << ", #bytes = " << nBytes << endl;

    char buffer[8];

    uint32_t* data = reinterpret_cast<uint32_t*>( buffer );

    data[0] = dataType;
    data[1] = nBytes;

    outFile.write( buffer, 8 );
}

/* --------------------------------------------------------------------------------- */

static inline uint32_t paddingBytes( const uint32_t nBytes )
{
    uint32_t padSize = nBytes % 8;

    if ( padSize )
    {
        padSize = 8 - padSize;
    }

    return padSize;
}

/* --------------------------------------------------------------------------------- */

static void writePadding( IOStream& outFile, uint32_t size )
{
    uint32_t padSize = paddingBytes( size );

    static char padData[] = { 0, 0, 0, 0, 0, 0, 0, 0 };

    outFile.write( padData, padSize );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType> 
static void writeData( IOStream& outFile, const ValueType* data, uint32_t size )
{
    uint32_t dataType  = matlabType<ValueType>();

    uint32_t nBytes    = size * sizeof( ValueType );

    writeDataElementHeader( outFile, dataType, nBytes );

    outFile.write( reinterpret_cast<const char*>( data ), nBytes );
    writePadding( outFile, nBytes );
}

/* --------------------------------------------------------------------------------- */

static void writeData( IOStream& outFile, const char* name )
{
    uint32_t size = strlen( name );

    writeData( outFile, name, size );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
static void writeData( IOStream& outFile, const HArray<ValueType>& array )
{
    if ( isComplex ( array.getValueType() ) )
    {
        typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

        HArray<AbsType> real;

        utilskernel::HArrayUtils::setArray( real, array );
        writeData( outFile, real );

        HArray<ValueType> tmp;  
        ValueType minusi = ComplexDouble( 0, -1 );
        utilskernel::HArrayUtils::binaryOpScalar2( tmp, array, minusi, utilskernel::binary::MULT );
        utilskernel::HArrayUtils::setArray( real, tmp );

        writeData( outFile, real );

        return;
    }

    std::cout << "writeData, array = " << array << ", size = " << array.size() << std::endl;

    uint32_t arraySize = array.size();

    {
        ReadAccess<ValueType> rArray( array );
        writeData( outFile, rArray.get(), arraySize );
    }
}

/* --------------------------------------------------------------------------------- */

static void writeSparseHeader( IOStream& outFile, const uint32_t n, const uint32_t m, const uint32_t nnz, const uint32_t nBytes, bool isComplex )
{
    uint32_t dataType = MAT_MATRIX;

    writeDataElementHeader( outFile, dataType, nBytes );

    uint32_t header[2];
    char*    headerBytes = reinterpret_cast<char*>( header );

    headerBytes[0] = MAT_SPARSE_CLASS;                // class
    headerBytes[1] = isComplex ? ( 1 << 3 ) : 0 ;     // array flags

    header[1] = nnz;                                  // # non-zero entries

    writeData( outFile, header, 2 );
  
    int dims[2] = { n, m };

    writeData( outFile, dims, 2 );
    writeData( outFile, "LAMA" );
}

/* --------------------------------------------------------------------------------- */

static void writeArrayHeader( IOStream& outFile, const IndexType n, const common::scalar::ScalarType stype )
{
    bool isComplex = common::isComplex( stype );
 
    uint32_t dataBytes = common::typeSize( stype );

    if ( isComplex ) 
    {
        dataBytes /= 2;   // real and imaginary part
    }

    dataBytes = n * dataBytes;

    uint32_t padSize = paddingBytes( dataBytes );

    if ( padSize )
    {
        dataBytes += padSize;
    }

    uint32_t dataType = MAT_MATRIX;
    uint32_t nBytes   = 16 + 16 + 16 + 8 + dataBytes;

    if ( isComplex )
    {
        nBytes += 8 + dataBytes;   // additional header for complex values
    }

    // Note: 16 bytes for flags, 16 bytes for dims, 16 bytes for name

    std::cout << "writeArrayHeader, n = " << n << ", type = " << stype << endl;

    writeDataElementHeader( outFile, dataType, nBytes );

    uint32_t header[2];
    char*    headerBytes = reinterpret_cast<char*>( header );

    memset( headerBytes, 0, 8 );

    uint8_t matClass = scalarType2Class( stype );

    headerBytes[0] = matClass;                        // class
    headerBytes[1] = isComplex ? ( 1 << 3 ) : 0 ;     // array flags

    std::cout << "Array flags, class = " << (int) headerBytes[0] << ", flags = " << (int) headerBytes[1] << std::endl;

    writeData( outFile, header, 2 );

    int dims[2];

    dims[0] = static_cast<int>( n );
    dims[1] = 1;

    writeData( outFile, dims, 2 );
    writeData( outFile, "LAMA" );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not supported for " << *this )

    IOStream outFile( fileName, std::ios::out );

    writeMATFileHeader( outFile );

    writeArrayHeader( outFile, array.size(), array.getValueType() );

    writeData( outFile, array );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not supported for " << *this )

    CSRStorage<ValueType> csrStorage;

    IndexType numRows = storage.getNumRows();
    IndexType numCols = storage.getNumColumns();

    csrStorage.assignTranspose( storage );

    IndexType numValues = csrStorage.getNumValues();

    IOStream outFile( fileName, std::ios::out );

    writeMATFileHeader( outFile );
    
    uint32_t size = 8;   //   number of bytes for element header

    size += 24;  // dims + name 
    size += 16 + numValues * sizeof( IndexType );
    size += 16 + ( numCols + 1 ) * sizeof( IndexType );
    size += 16 + numValues * sizeof( ValueType );

    writeSparseHeader( outFile, numRows, numCols, numValues, size, false );

    writeData( outFile, csrStorage.getJA() );
    writeData( outFile, csrStorage.getIA() );
    writeData( outFile, csrStorage.getValues() );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readData(
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

void MatlabIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName )
{
    IOStream inFile( fileName, std::ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;
    readMATFileHeader( version, endian, inFile );

    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    common::scoped_array<char> dataElement;

    readDataElement( inFile, dataElement );

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, dataElement.get() );

    SCAI_ASSERT_EQ_ERROR( dataType, MAT_MATRIX, "can only read storage as MATRIX - MATLAB array" )

    uint32_t header[2];
    int dims[2];

    uint32_t offset = getData( header, 2, dataPtr );
    offset += getData( dims, 2, dataPtr + offset );

    numValues  = header[1];
    numRows    = dims[0];
    numColumns = dims[1];
}

/* --------------------------------------------------------------------------------- */

static void getStorage( _MatrixStorage& storage, const char* dataElementPtr, uint32_t nBytes )
{
    std::cout << "getStorage, nBytes = " << nBytes << std::endl;

    // read subelements of miMatrix

    uint32_t offset = 0;

    uint32_t subType;
    uint32_t subNBytes;
    uint32_t subWBytes;

    readDataElementHeader( subType, subNBytes, subWBytes, dataElementPtr + offset );

    SCAI_ASSERT_EQ_ERROR( subNBytes, 8, "ArrayFlags must be 8 bytes" )
    SCAI_ASSERT_EQ_ERROR( subType, 6, "ArrayFlags type must be miUINT32" )

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

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readStorageImpl(
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

    readDataElement( inFile, dataElement );

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, dataElement.get() );
 
    SCAI_ASSERT_EQ_ERROR( dataType, MAT_MATRIX, "can only read storage as MATRIX - MATLAB array" )

    getStorage( storage, dataPtr, nBytes );

    SCAI_LOG_ERROR( logger, "readStorage: " << storage << ", #bytes read = " << wBytes << " / " << nBytes )
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
