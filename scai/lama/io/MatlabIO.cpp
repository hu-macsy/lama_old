/**
 * @file MatlabIO.cpp
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
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/io/MATIOStream.hpp>
#include <scai/lama/io/ImageIO.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/exception/IOException.hpp>

#include <sstream>

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

FileIO* MatlabIO::create()
{
    return new MatlabIO();
}

string MatlabIO::createValue()
{
    return MAT_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool MatlabIO::isSupportedMode( const FileMode mode ) const
{
    // only binary is supported

    if ( mode == FORMATTED )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeAt( ostream& stream ) const
{
    stream << "MatlabIO ( suffix = " << MAT_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatlabIO::logger, "FileIO.MatlabIO" )

/* --------------------------------------------------------------------------------- */

template<typename ArrayType, typename DataType>
void MatlabIO::readMATArrayImpl( hmemo::HArray<ArrayType>& array, const void* data, IndexType nBytes )
{
    IndexType elemSize  = sizeof( DataType );
    IndexType arraySize = nBytes / elemSize;

    SCAI_ASSERT_EQ_ERROR( elemSize* arraySize, nBytes, "Size mismatch, elemSize = " << elemSize << ", arraySize = " << arraySize )

    if ( typeid( ArrayType ) == typeid( DataType ) )
    {
        SCAI_LOG_INFO( logger, "readMATArrayImpl, in place, type = " << common::TypeTraits<ArrayType>::id()
                       << ", arraySize = " << arraySize << ", #bytes = " << nBytes )

        // no temporary array required

        hmemo::WriteOnlyAccess<ArrayType> wData( array, arraySize );
        ::memcpy( wData.get(), data, nBytes );

    }
    else
    {
        SCAI_LOG_INFO( logger, "readMATArrayImpl, in place, type = " << common::TypeTraits<ArrayType>::id()
                       << ", arraySize = " << arraySize << ", #bytes = " << nBytes )

        // temporary array and conversion required

        hmemo::HArray<DataType> tmp;

        {
            hmemo::WriteOnlyAccess<DataType> wData( tmp, arraySize );
            ::memcpy( wData.get(), data, nBytes );
        }

        HArrayUtils::assign( array, tmp );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readMATArray( hmemo::HArray<ValueType>& array, const char* data, const uint32_t mxType, const uint32_t nBytes )
{
    switch ( mxType )
    {
        case MATIOStream::MAT_DOUBLE  :
            readMATArrayImpl<ValueType, double>( array, data, nBytes );
            break;
        case MATIOStream::MAT_FLOAT   :
            readMATArrayImpl<ValueType, float>( array, data, nBytes );
            break;
        case MATIOStream::MAT_LDOUBLE :
            readMATArrayImpl<ValueType, long double>( array, data, nBytes );
            break;
        case MATIOStream::MAT_INT8    :
            readMATArrayImpl<ValueType, int8_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_UINT8   :
            readMATArrayImpl<ValueType, uint8_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_INT16   :
            readMATArrayImpl<ValueType, int16_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_UINT16  :
            readMATArrayImpl<ValueType, uint16_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_INT32   :
            readMATArrayImpl<ValueType, int32_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_UINT32  :
            readMATArrayImpl<ValueType, uint32_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_INT64   :
            readMATArrayImpl<ValueType, int64_t>( array, data, nBytes );
            break;
        case MATIOStream::MAT_UINT64  :
            readMATArrayImpl<ValueType, uint64_t>( array, data, nBytes );
            break;

        default :
            COMMON_THROWEXCEPTION( "mxType = " << mxType << " is unknown data type in Matlab file." )
    }
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readArrayInfo( IndexType& n, const string& arrayFileName )
{
    MATIOStream inFile( arrayFileName, ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    inFile.readMATFileHeader( version, endian );

    common::scoped_array<char> dataElement;

    inFile.readDataElement( dataElement );

    IndexType nDims;
    IndexType dims[2];
    IndexType nnz;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    MATIOStream::getMatrixInfo( matClass, dims, 2, nDims, nnz, isComplex, dataElement.get() );

    n = dims[0] * dims[1];

    if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
    {
        COMMON_THROWEXCEPTION( "File " << arrayFileName << " contains sparse matrix, but not array" )
    }

    if ( MATIOStream::class2ScalarType( matClass ) == common::scalar::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File " << arrayFileName << " contains unsupported matrix class = " << matClass )
    }

    if ( dims[1] != 1 || dims[0] != 1 )
    {
        SCAI_LOG_WARN( logger, "File " << arrayFileName << ": matrix " << dims[0] << " x " << dims[1] << " considered as array of size " << n )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void buildComplex( HArray<ValueType>& array, HArray<ValueType>& imagValues )
{
    if ( !common::isComplex( array.getValueType() ) )
    {
        // SCAI_LOG_WARN( logger, "imaginary values are ignored" )
        return;
    }

#ifdef SCAI_COMPLEX_SUPPORTED

    // array = array + i * imagValues

    ValueType i = static_cast<ValueType>( ComplexDouble( 0, 1 ) );
    utilskernel::HArrayUtils::compute( imagValues, imagValues, common::binary::MULT, i );
    utilskernel::HArrayUtils::binaryOp( array, array, imagValues, common::binary::ADD );

#endif
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const string& arrayFileName,
    const IndexType ,
    const IndexType )
{
    MATIOStream inFile( arrayFileName, ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    inFile.readMATFileHeader( version, endian );

    common::scoped_array<char> dataElement;

    uint32_t nBytes = inFile.readDataElement( dataElement );

    SCAI_LOG_INFO( logger, "Read full data element from input file, #bytes = " << nBytes )

    const char* elementPtr = dataElement.get();

    IndexType dims[2];
    IndexType nnz;
    IndexType nDims;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    uint32_t offset = MATIOStream::getMatrixInfo( matClass, dims, 2, nDims, nnz, isComplex, elementPtr );

    SCAI_ASSERT_LE_ERROR( offset, nBytes, "data element insufficient to read matrix info" )

    // now read the data

    offset += getArrayData( array, elementPtr + offset, nBytes - offset );

    SCAI_ASSERT_EQ_ERROR( array.size(), dims[0] * dims[1], "serious mismatch" )

    if ( isComplex )
    {
        utilskernel::LArray<ValueType> imagValues;

        offset += getArrayData( imagValues, elementPtr + offset, nBytes - offset );

        buildComplex( array, imagValues );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const std::string& fileName )
{
    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray, fileName, 0, nIndex );
    size = denseArray.size();
    utilskernel::HArrayUtils::buildSparseArrayImpl( values, indexes, denseArray );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
uint32_t MatlabIO::writeArrayData( MATIOStream& outFile, const HArray<ValueType>& array, bool dryRun )
{
    uint32_t wBytes = 0;

    if ( isComplex ( array.getValueType() ) )
    {
        typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

        HArray<AbsType> real;

        utilskernel::HArrayUtils::setArray( real, array );

        wBytes += writeArrayData( outFile, real, dryRun );

        HArray<ValueType> tmp;
        ValueType minusi = ComplexDouble( 0, -1 );
        utilskernel::HArrayUtils::compute( tmp, array, common::binary::MULT, minusi );
        utilskernel::HArrayUtils::setArray( real, tmp );

        wBytes += writeArrayData( outFile, real, dryRun );
    }
    else
    {
        SCAI_LOG_INFO( logger, "writeArrayData, array = " << array << ", size = " << array.size() )

        uint32_t arraySize = array.size();

        {
            ReadAccess<ValueType> rArray( array );
            wBytes = outFile.writeData( rArray.get(), arraySize, dryRun );
        }
    }

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeDenseArray( MATIOStream& outFile, const hmemo::HArray<ValueType>& array, IndexType dims[] )
{
    SCAI_ASSERT_EQ_ERROR( array.size(), dims[0] * dims[1], "array size / dims mismatch" )

    common::scalar::ScalarType stype = array.getValueType();

    if ( mScalarTypeData == common::scalar::PATTERN )
    {
        COMMON_THROWEXCEPTION( "Cannot write dense data as pattern" )
    }

    if ( mScalarTypeData != common::scalar::INTERNAL )
    {
        if ( stype != mScalarTypeData )
        {
            SCAI_LOG_WARN( logger, "write HArray<" << stype << "> as it is, IO_TYPE=" << mScalarTypeData << " ignored" )
        }
    }

    uint32_t nBytes = 16;  // initial guess used for the dry run

    bool dryRun = true;   // make a dryRun at first to determine the size of written bytes

    uint32_t wBytes = outFile.writeDenseHeader( dims[0], dims[1], nBytes, stype, dryRun );
    wBytes += writeArrayData( outFile, array, dryRun );

    nBytes = wBytes - 8;  // subtract for the first header

    SCAI_LOG_INFO( logger, "writeDenseArray, dryrun gives written bytes = " << wBytes << ", now write" )

    dryRun = false;  // now write it with the correct value of nBytes

    wBytes  = outFile.writeDenseHeader( dims[0], dims[1], nBytes, stype, dryRun );
    wBytes += writeArrayData( outFile, array, dryRun );

    SCAI_LOG_INFO( logger, "written dense array " << array << " as " << dims[0] << " x " << dims[1]
                   << ", wBytes = " << wBytes )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const string& fileName )
{
    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not supported for " << *this )

    MATIOStream outFile( fileName, ios::out );

    outFile.writeMATFileHeader();

    IndexType dims[2] = { array.size(), 1 };

    writeDenseArray( outFile, array, dims );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeSparseImpl(
    const IndexType size,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values,
    const std::string& fileName )
{
    // sparse unsupported for this file format, write it dense

    HArray<ValueType> denseArray;
    utilskernel::HArrayUtils::buildDenseArray( denseArray, size, values, indexes );
    writeArrayImpl( denseArray, fileName );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const string& fileName )
{
    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not supported for " << *this )

    IndexType numRows = storage.getNumRows();
    IndexType numCols = storage.getNumColumns();
    IndexType numValues = storage.getNumValues();

    MATIOStream outFile( fileName, ios::out );

    outFile.writeMATFileHeader();

    if ( numValues * 2 >= numRows* numCols && mScalarTypeData != common::scalar::PATTERN )
    {
        SCAI_LOG_INFO( logger, "Write storage as dense matrix to file " << fileName << ": " << storage )

        DenseStorage<ValueType> denseStorage;

        denseStorage.assignTranspose( storage );   // MATLAB stores it column-wise

        HArray<ValueType>& array = denseStorage.getData();

        IndexType dims[2] = { numRows, numCols };

        writeDenseArray( outFile, array, dims );
    }
    else
    {
        SCAI_LOG_INFO( logger, "Write storage as sparse matrix to file " << fileName << ": " << storage )

        HArray<IndexType> ia;
        HArray<IndexType> ja;
        HArray<ValueType> values;

        storage.buildCSCData( ja, ia, values );

        IndexType numValues = storage.getNumValues();

        uint32_t wBytes = 24;   // initial guess for dryrun, avoid short write

        bool isComplex = common::isComplex( storage.getValueType() );

        for ( int i = 0; i < 2; ++i )
        {
            bool dryRun = ( i == 0 );      // first run dry, second run okay

            wBytes  = outFile.writeSparseHeader( numRows, numCols, numValues, wBytes, isComplex, dryRun );
            wBytes += writeArrayData( outFile, ia, dryRun );
            wBytes += writeArrayData( outFile, ja, dryRun );

            if ( mScalarTypeData != common::scalar::PATTERN )
            {
                wBytes += writeArrayData( outFile, values, dryRun );
            }

            wBytes -= 8;  // subtract for the first header

            SCAI_LOG_INFO( logger, "writeStorage, dryRun = " << dryRun << ", wBytes = " << wBytes )
        }
    }
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const string& fileName )
{
    MATIOStream inFile( fileName, ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    inFile.readMATFileHeader( version, endian );

    common::scoped_array<char> dataElement;

    inFile.readDataElement( dataElement );

    IndexType dims[2];
    IndexType nDims;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    MATIOStream::getMatrixInfo( matClass, dims, 2, nDims, numValues, isComplex, dataElement.get() );

    numRows    = dims[0];
    numColumns = dims[1];

    if ( matClass != MATIOStream::MAT_SPARSE_CLASS )
    {
        if ( MATIOStream::class2ScalarType( matClass ) == common::scalar::UNKNOWN )
        {
            COMMON_THROWEXCEPTION( "File " << fileName << " contains unsupported matrix class = " << matClass )
        }

        numValues  = dims[0] * dims[1];
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
uint32_t MatlabIO::getArrayData( HArray<ValueType>& array, const char* data, uint32_t len )
{
    SCAI_ASSERT_LE_ERROR( 8, len, "data insufficient to read element header" )

    uint32_t wBytes;
    uint32_t nBytes;
    uint32_t dataType;

    const char* arrayDataPtr = MATIOStream::readDataElementHeader( dataType, nBytes, wBytes, data );

    SCAI_ASSERT_LE_ERROR( wBytes, len, "data insufficient to read array" )

    readMATArray( array, arrayDataPtr, dataType, nBytes );

    SCAI_LOG_INFO( logger, "read array " << array << " from data ( len = " << len << " ), type = " << dataType
                   << " is " << MATIOStream::matlabType2ScalarType( dataType )
                   << ", nBytes = " << nBytes << ", wBytes = " << wBytes )

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

template <typename ValueType>
uint32_t MatlabIO::getSparseStorage( MatrixStorage<ValueType>& storage,
                                     const IndexType dims[2], const IndexType nnz,
                                     bool isComplex,
                                     const char* dataElementPtr, uint32_t nBytes )
{
    SCAI_LOG_INFO( logger, "Get sparse<" << common::TypeTraits<ValueType>::stype << "> matrix "
                   << dims[0] << " x " << dims[1] << ", nnz = " << nnz )

    uint32_t offset = 0;

    // read IA, JA, Values of sparse array

    HArray<IndexType> ia;
    HArray<IndexType> ja;
    HArray<ValueType> values;

    offset += getArrayData( ia, dataElementPtr + offset, nBytes - offset );
    offset += getArrayData( ja, dataElementPtr + offset, nBytes - offset );

    if ( mScalarTypeData == common::scalar::PATTERN )
    {
        values.init( ValueType( 1 ), nnz );   // set values with default value
    }
    else
    {
        offset += getArrayData( values, dataElementPtr + offset, nBytes - offset );

        if ( isComplex )
        {
            HArray<ValueType> imagValues;
            offset += getArrayData( imagValues, dataElementPtr + offset, nBytes - offset );
            buildComplex( values, imagValues );
        }
    }

    CSRStorage<ValueType> csrStorage;
    csrStorage.allocate( dims[1], dims[0] );  // will be transposed
    csrStorage.swap( ja, ia, values );
    csrStorage.assignTranspose( csrStorage );
    csrStorage.sortRows( dims[0] == dims[1] );
    storage = csrStorage;

    return offset;
}

/* --------------------------------------------------------------------------------- */

template <typename ValueType>
uint32_t MatlabIO::getStructStorage( MatrixStorage<ValueType>& storage, const char* dataElementPtr, uint32_t nBytes )
{
    int len;
    char names[1024];

    uint32_t offset = MATIOStream::getData( &len, 1, dataElementPtr );
    uint32_t size   = MATIOStream::getString( names, 1024, dataElementPtr + offset );

    offset += size;

    int nFields = ( size - 8 ) / len;

    SCAI_LOG_INFO( logger, "parse structure with " << nFields << " fields, name len = " << len )

    bool readStorage = false;    // set it to true if any sparse array is found

    for ( int i = 0; i < nFields; ++i )
    {
        SCAI_ASSERT_LT_ERROR( offset, nBytes, "no more data for fields" )

        char* ptr = names + i * len;

        SCAI_LOG_INFO( logger, "parse field " << i << " of " << nFields << ": name = " << i << " = " << ptr )

        uint32_t dataType;
        uint32_t nBytesField;
        uint32_t wBytes;

        MATIOStream::readDataElementHeader( dataType, nBytesField, wBytes, dataElementPtr + offset );

        IndexType dims[2];

        IndexType nnz;
        IndexType nDims;
        IndexType maxDims = 2;
        bool      isComplex;
        MATIOStream::MATClass matClass;

        SCAI_LOG_INFO( logger, "read structure field[" << i << "], name = " << ptr << ", dataType = " << dataType
                       << ", nBytes = " << nBytesField << " / " << wBytes )

        uint32_t offset1 = 0;

        nBytesField = wBytes; // reset it

        offset1 += MATIOStream::getMatrixInfo( matClass, dims, maxDims, nDims, nnz, isComplex, dataElementPtr + offset + offset1, false );

        SCAI_LOG_INFO( logger, "read info of cell " << dims[0] << " x " << dims[1]
                       << ", nnz = " << nnz << ", isComplex = " << isComplex << ", class = " << matClass )

        if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
        {
            SCAI_ASSERT_ERROR( !readStorage, "more than one sparse array identified in struct" )
            getSparseStorage( storage, dims, nnz, isComplex, dataElementPtr + offset + offset1, nBytesField - offset1 );
            SCAI_LOG_INFO( logger, "read sparse storage = " << storage )
            readStorage = true;
            // continue loop to parse all further elements
        }

        offset += wBytes;
    }

    SCAI_ASSERT_ERROR( readStorage, "no storage found in struct with " << nFields << " fields" )

    return offset;
}

/* --------------------------------------------------------------------------------- */

template <typename ValueType>
void MatlabIO::getStorage( MatrixStorage<ValueType>& storage, const char* dataElementPtr, uint32_t nBytes )
{
    IndexType dims[2];
    IndexType nnz;
    IndexType nDims;
    const IndexType maxDims = 2;
    bool      isComplex;
    MATIOStream::MATClass matClass;

    uint32_t offset = MATIOStream::getMatrixInfo( matClass, dims, maxDims, nDims, nnz, isComplex, dataElementPtr );

    if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
    {
        offset += getSparseStorage( storage, dims, nnz, isComplex, dataElementPtr + offset, nBytes - offset );
    }
    else if ( matClass == MATIOStream::MAT_STRUCT_CLASS )
    {
        // dims = [1, 1]
        offset += getStructStorage( storage, dataElementPtr + offset, nBytes - offset );
    }
    else if ( matClass == MATIOStream::MAT_CELL_CLASS )
    {
        COMMON_THROWEXCEPTION( "Cell Array Data Element Format not supported yet" )
    }
    else if ( matClass == MATIOStream::MAT_OBJECT_CLASS )
    {
        COMMON_THROWEXCEPTION( "Object MAT-File Data Element Format not supported yet" )
    }
    else if ( MATIOStream::class2ScalarType( matClass ) == common::scalar::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File contains unsupported matrix class = " << matClass )
    }
    else
    {
        SCAI_LOG_INFO( logger, "Get dense<" << common::TypeTraits<ValueType>::stype << "> matrix "
                       << dims[0] << " x " << dims[1] )

        LArray<ValueType> values;

        offset += getArrayData( values, dataElementPtr + offset, nBytes - offset );

        if ( isComplex )
        {
            HArray<ValueType> imagValues;
            offset += getArrayData( imagValues, dataElementPtr + offset, nBytes - offset );
            buildComplex( values, imagValues );
        }

        // MATLAB stores it columnwise, so we transpose the data

        storage.setDenseData( dims[1], dims[0], values );
        storage.assignTranspose( storage );
    }

    SCAI_ASSERT_EQ_ERROR( offset, nBytes, "mismatch read bytes and size bytes, maybe COMPLEX" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const string& matrixFileName,
    const IndexType firstRow,
    const IndexType nRows )
{
    MATIOStream inFile( matrixFileName, ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    inFile.readMATFileHeader( version, endian );

    common::scoped_array<char> dataElement;

    uint32_t nBytes = inFile.readDataElement( dataElement );

    if ( firstRow == 0 && nRows == nIndex )
    {
        getStorage( storage, dataElement.get(), nBytes );
        SCAI_LOG_INFO( logger, "readStorage: " << storage )
    }
    else
    {
        COOStorage<ValueType> tmpStorage;
        getStorage( tmpStorage, dataElement.get(), nBytes );
        SCAI_LOG_INFO( logger, "readStorage: " << tmpStorage )
        tmpStorage.copyBlockTo( storage, firstRow, nRows );
        SCAI_LOG_INFO( logger, "extracted: first = " << firstRow << ", #rows = " << nRows << ": " << storage )
    }
}

template<typename ValueType>
static void changeMajor( hmemo::HArray<ValueType>& out, const hmemo::HArray<ValueType>& in, const common::Grid& grid )
{
    SCAI_ASSERT_EQ_ERROR( 3, grid.nDims(), "other dims not supported yet" )
    hmemo::ReadAccess<ValueType> rIn( in );
    hmemo::WriteOnlyAccess<ValueType> wOut( out, grid.size() );
    const IndexType n0 = grid.size(0);
    const IndexType n1 = grid.size(1);
    const IndexType n2 = grid.size(2);

    const IndexType dIn0 = 1;
    const IndexType dIn1 = n0;
    const IndexType dIn2 = n0 * n1;
    const IndexType dOut0 = n1 * n2;
    const IndexType dOut1 = n2;
    const IndexType dOut2 = 1;

    for ( IndexType i0 = 0; i0 < n0; ++i0 )
    for ( IndexType i1 = 0; i1 < n1; ++i1 )
    for ( IndexType i2 = 0; i2 < n2; ++i2 )
    {
        wOut[ i0 * dOut0 + i1 * dOut1 + i2 * dOut2 ] = rIn[ i0 * dIn0 + i1 * dIn1 + i2 * dIn2 ];
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Read grid vector data from file                                              */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void MatlabIO::readImpl( HArray<ValueType>& data, common::Grid& grid, const std::string& gridFileName )
{
    MATIOStream inFile( gridFileName, ios::in );

    int version = 0;
    IOStream::Endian endian = IOStream::MACHINE_ENDIAN;

    inFile.readMATFileHeader( version, endian );

    common::scoped_array<char> dataElement;

    uint32_t nBytes = inFile.readDataElement( dataElement );

    const char* elementPtr = dataElement.get();

    IndexType nDims;
    IndexType dims[SCAI_GRID_MAX_DIMENSION];
    IndexType nnz;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    uint32_t offset = MATIOStream::getMatrixInfo( matClass, dims, SCAI_GRID_MAX_DIMENSION, nDims, nnz, isComplex, dataElement.get() );

    if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
    {
        COMMON_THROWEXCEPTION( "File " << gridFileName << " contains sparse matrix, but not grid array" )
    }

    if ( MATIOStream::class2ScalarType( matClass ) == common::scalar::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File " << gridFileName << " contains unsupported matrix class = " << matClass )
    }

    grid = common::Grid( nDims, dims );

    SCAI_ASSERT_LE_ERROR( offset, nBytes, "data element insufficient to read matrix info" )

    // now read the data

    offset += getArrayData( data, elementPtr + offset, nBytes - offset );

    SCAI_ASSERT_EQ_ERROR( data.size(), grid.size(), "serious mismatch" )

    if ( isComplex )
    {
        utilskernel::LArray<ValueType> imagValues;

        offset += getArrayData( imagValues, elementPtr + offset, nBytes - offset );

        buildComplex( data, imagValues );
    }

    HArray<ValueType> tmpData( grid.size() );
    changeMajor( tmpData, data, grid );
    data.swap( tmpData );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readGridArray( _HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::read( ( MatlabIO& ) *this, data, grid, inputFileName );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    SCAI_ASSERT_EQ_ERROR( data.size(), grid.size(), "size of array does not match the grid size" )

    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to Matlab file" )
    }

    CRTPFileIO<MatlabIO>::writeArray( data, outputFileName );
}

/* --------------------------------------------------------------------------------- */

#define SCAI_MATLAB_METHOD_INSTANTIATIONS( _type )           \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::writeArrayImpl(                           \
        const hmemo::HArray<_type>& array,                   \
        const string& fileName );                            \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::readArrayImpl(                            \
        hmemo::HArray<_type>& array,                         \
        const string& arrayFileName,                         \
        const IndexType ,                                    \
        const IndexType );                                   \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::writeSparseImpl(                          \
        const IndexType size,                                \
        const HArray<IndexType>& index,                      \
        const HArray<_type>& values,                         \
        const std::string& fileName );                       \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::readSparseImpl(                           \
        IndexType& size,                                     \
        HArray<IndexType>& indexes,                          \
        HArray<_type>& values,                               \
        const std::string& fileName );         

SCAI_COMMON_LOOP( SCAI_MATLAB_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_MATLAB_METHOD_INSTANTIATIONS

#define SCAI_MATLAB_METHOD_INSTANTIATIONS( _type )      \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void MatlabIO::writeStorageImpl(                    \
        const MatrixStorage<_type>& storage,            \
        const string& fileName );                       \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void MatlabIO::readStorageImpl(                     \
        MatrixStorage<_type>& storage,                  \
        const string& matrixFileName,                   \
        const IndexType firstRow,                       \
        const IndexType nRows );                     

SCAI_COMMON_LOOP( SCAI_MATLAB_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_MATLAB_METHOD_INSTANTIATIONS

}  // lama

}  // scai


