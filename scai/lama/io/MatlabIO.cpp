/**
 * @file MatlabIO.cpp
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
 * @brief Implementation of methods for FileIO class MatlabIO
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include <scai/lama/io/MatlabIO.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/io/MATIOStream.hpp>
#include <scai/lama/io/ImageIO.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/utilskernel/openmp/OpenMPSection.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/exception/IOException.hpp>
#include <scai/common/safer_memcpy.hpp>

#include <sstream>
#include <memory>

#define MAT_SUFFIX ".mat"

using scai::common::safer_memcpy;

using namespace std;

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatlabIO::logger, "FileIO.MatlabIO" )

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

    if ( mode == FileMode::FORMATTED )
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

void MatlabIO::open( const char* fileName, const char* fileMode )
{
    if ( strcmp( fileMode, "w" ) == 0 )
    {
        SCAI_ASSERT( mFileMode != FileMode::FORMATTED, "Formatted output not supported for " << *this )

        mFile.open( fileName, std::ios::out | std::ios::trunc );
        mFile.writeMATFileHeader();
    }
    else if ( strcmp( fileMode, "r" ) == 0 )
    {
        mFile.open( fileName, std::ios::in );
        int version = 0;
        IOStream::Endian endian = IOStream::MACHINE_ENDIAN;
        mFile.readMATFileHeader( version, endian );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported file mode for Matlab file: " << fileMode )
    }
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::close()
{
    mFile.close();
}

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
        safer_memcpy( wData.get(), data, nBytes );

    }
    else
    {
        SCAI_LOG_INFO( logger, "readMATArrayImpl, in place, type = " << common::TypeTraits<ArrayType>::id()
                       << ", arraySize = " << arraySize << ", #bytes = " << nBytes )

        // temporary array and conversion required

        hmemo::HArray<DataType> tmp;

        {
            hmemo::WriteOnlyAccess<DataType> wData( tmp, arraySize );
            safer_memcpy( wData.get(), data, nBytes );
        }

        // Attention: not all conversions might be supported

        HArrayUtils::_assign( array, tmp );
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

void MatlabIO::getArrayInfo( IndexType& n )
{
    std::unique_ptr<char[]> dataElement;

    std::streampos pos = mFile.tellg();

    mFile.readDataElement( dataElement );

    mFile.clear();       // important to reset flags
    mFile.seekg( pos );

    IndexType nDims;
    IndexType dims[2];
    IndexType nnz;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    MATIOStream::getMatrixInfo( matClass, dims, 2, nDims, nnz, isComplex, dataElement.get() );

    n = dims[0] * dims[1];

    if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
    {
        COMMON_THROWEXCEPTION( "File " << mFile.getFileName() << " contains sparse matrix, but not array" )
    }

    if ( MATIOStream::class2ScalarType( matClass ) == common::ScalarType::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File " << mFile.getFileName() << " contains unsupported matrix class = " << matClass )
    }

    if ( dims[1] != 1 || dims[0] != 1 )
    {
        SCAI_LOG_WARN( logger, "File " << mFile.getFileName() << ": matrix " << dims[0] << " x " << dims[1] << " considered as array of size " << n )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void buildComplex( HArray<ValueType>& array, HArray<ValueType>& imagValues )
{
    if ( !common::isComplex( array.getValueType() ) )
    {
        return;
    }


    // array = array + i * imagValues

#ifdef SCAI_COMPLEX_SUPPORTED
    ValueType i = static_cast<ValueType>( ComplexDouble( 0, 1 ) );
#else
    ValueType i = -1;    // Note: will never be called here
#endif

    utilskernel::HArrayUtils::compute( imagValues, imagValues, common::BinaryOp::MULT, i );
    utilskernel::HArrayUtils::binaryOp( array, array, imagValues, common::BinaryOp::ADD );

}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
static void changeMajor( hmemo::HArray<ValueType>& out, const hmemo::HArray<ValueType>& in, const common::Grid& grid, bool isRow2Col )
{
    IndexType nDims = grid.nDims();

    // arrays for distances between elements, one for row-major ordering, one for col-major

    IndexType rowMajorDist[ SCAI_GRID_MAX_DIMENSION ];  // used in LAMA
    IndexType colMajorDist[ SCAI_GRID_MAX_DIMENSION ];  // required for MATLAB

    grid.getDistances( rowMajorDist );  

    // column major ordering of grid( n1, n2, n3, n4 ) is dist (1, n1, n1*n2, n1*n2*n3 )

    if ( nDims > SCAI_GRID_MAX_DIMENSION )
    {
        COMMON_THROWEXCEPTION( "serious: too many dims" )
    }
    else
    {
        colMajorDist[ 0 ] = 1;

        for ( IndexType i = 1; i < nDims; ++i )
        {
            colMajorDist[i] = colMajorDist[ i - 1 ] * grid.size( i - 1 );
        }
    }

    hmemo::ReadAccess<ValueType> rIn( in );
    hmemo::WriteOnlyAccess<ValueType> wOut( out, grid.size() );

    if ( isRow2Col )
    {
        // convert row-major ordering to column-major ordering

        utilskernel::OpenMPSection::assign( wOut.get(), nDims, grid.sizes(), colMajorDist, 
                                            rIn.get(), rowMajorDist, 
                                            common::BinaryOp::COPY, false );
    }
    else
    {
        // convert column-major ordering to row-major ordering

        utilskernel::OpenMPSection::assign( wOut.get(), nDims, grid.sizes(), rowMajorDist, 
                                            rIn.get(), colMajorDist, 
                                            common::BinaryOp::COPY, false );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readArrayImpl( hmemo::HArray<ValueType>& array )
{
    std::unique_ptr<char[]> dataElement;

    uint32_t nBytes = mFile.readDataElement( dataElement );

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
        HArray<ValueType> imagValues;

        offset += getArrayData( imagValues, elementPtr + offset, nBytes - offset );

        buildComplex( array, imagValues );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values )
{
    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray );
    size = denseArray.size();
    ValueType zeroValue = 0;
    utilskernel::HArrayUtils::buildSparseArray( values, indexes, denseArray, zeroValue );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
uint32_t MatlabIO::writeArrayData( const HArray<ValueType>& array, bool dryRun )
{
    uint32_t wBytes = 0;

    if ( isComplex ( array.getValueType() ) )
    {
        typedef typename common::TypeTraits<ValueType>::RealType RealType;

        HArray<RealType> real;

        utilskernel::HArrayUtils::setArray( real, array );

        wBytes += writeArrayData( real, dryRun );

        HArray<ValueType> tmp;
#ifdef SCAI_COMPLEX_SUPPORTED
        ValueType minusi = ComplexDouble( 0, -1 );
#else
        ValueType minusi = -1;  // unreachable code at all here
#endif
        utilskernel::HArrayUtils::compute( tmp, array, common::BinaryOp::MULT, minusi );
        utilskernel::HArrayUtils::setArray( real, tmp );

        wBytes += writeArrayData( real, dryRun );
    }
    else
    {
        SCAI_LOG_INFO( logger, "writeArrayData, array = " << array << ", size = " << array.size() )

        uint32_t arraySize = array.size();

        {
            ReadAccess<ValueType> rArray( array );
            wBytes = mFile.writeData( rArray.get(), arraySize, dryRun );
        }
    }

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeDenseGrid( const hmemo::HArray<ValueType>& array, const common::Grid& grid )
{
    SCAI_ASSERT_EQ_ERROR( array.size(), grid.size(), "array size / dims mismatch" )

    common::ScalarType stype = array.getValueType();

    if ( mScalarTypeData == common::ScalarType::PATTERN )
    {
        COMMON_THROWEXCEPTION( "Cannot write data as pattern" )
    }

    if ( mScalarTypeData != common::ScalarType::INTERNAL )
    {
        if ( stype != mScalarTypeData )
        {
            SCAI_LOG_WARN( logger, "write HArray<" << stype << "> as it is, IO_TYPE=" << mScalarTypeData << " ignored" )
        }
    }

    uint32_t nBytes = 16;  // initial guess used for the dry run

    bool dryRun = true;   // make a dryRun at first to determine the size of written bytes

    uint32_t wBytes = mFile.writeShapeHeader( grid.sizes(), grid.nDims(), nBytes, stype, dryRun );

    wBytes += writeArrayData( array, dryRun );

    nBytes = wBytes - 8;  // subtract for the first header

    SCAI_LOG_INFO( logger, "writeDenseGrid, dryrun gives written bytes = " << wBytes << ", now write" )

    dryRun = false;  // now write it with the correct value of nBytes

    wBytes  = mFile.writeShapeHeader( grid.sizes(), grid.nDims(), nBytes, stype, dryRun );
    wBytes += writeArrayData( array, dryRun );

    SCAI_LOG_INFO( logger, "written shaped array " << array << " with shape " << grid
                   << ", wBytes = " << wBytes )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeArrayImpl( const hmemo::HArray<ValueType>& array )
{
    // Matlab expects even for a vector a two-dimensional grid.

    const common::Grid2D grid( array.size(), 1 );

    writeDenseGrid( array, grid );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid )
{
    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeGridImpl( *this, data, grid );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeGridImpl(
    const hmemo::HArray<ValueType>& data,
    const common::Grid& grid )
{
    SCAI_LOG_INFO( logger, "writeGridImpl<" << common::TypeTraits<ValueType>::id() 
                            << ">, shape = " << grid << ", data = " << data )

    HArray<ValueType> tmpData( grid.size() );
 
    bool isRow2Col = true;

    changeMajor( tmpData, data, grid, isRow2Col );

    writeDenseGrid( tmpData, grid );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeSparseImpl(
    const IndexType size,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values )
{
    // sparse unsupported for this file format, write it dense

    HArray<ValueType> denseArray;
    ValueType zero = 0;
    utilskernel::HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );
    writeArrayImpl( denseArray );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeStorageImpl( const MatrixStorage<ValueType>& storage )
{
    IndexType numRows = storage.getNumRows();
    IndexType numCols = storage.getNumColumns();
    IndexType numValues = storage.getNumValues();

    if ( numValues * 2 >= numRows* numCols && mScalarTypeData != common::ScalarType::PATTERN )
    {
        SCAI_LOG_INFO( logger, "Write storage as dense matrix to file " << mFile.getFileName() << ": " << storage )

        DenseStorage<ValueType> denseStorage;

        denseStorage.assignTranspose( storage );   // MATLAB stores it column-wise

        const HArray<ValueType>& array = denseStorage.getValues();

        common::Grid2D grid( numRows, numCols );

        writeDenseGrid( array, grid );
    }
    else
    {
        SCAI_LOG_INFO( logger, "Write storage as sparse matrix to file " << mFile.getFileName() << ": " << storage )

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

            wBytes  = mFile.writeSparseHeader( numRows, numCols, numValues, wBytes, isComplex, dryRun );
            wBytes += writeArrayData( ia, dryRun );
            wBytes += writeArrayData( ja, dryRun );

            if ( mScalarTypeData != common::ScalarType::PATTERN )
            {
                wBytes += writeArrayData( values, dryRun );
            }

            wBytes -= 8;  // subtract for the first header

            SCAI_LOG_INFO( logger, "writeStorage, dryRun = " << dryRun << ", wBytes = " << wBytes )
        }
    }
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    std::unique_ptr<char[]> dataElement;

    std::streampos pos = mFile.tellg();

    mFile.readDataElement( dataElement );

    mFile.clear();       // important to reset flags
    mFile.seekg( pos );

    IndexType dims[2];
    IndexType nDims;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    MATIOStream::getMatrixInfo( matClass, dims, 2, nDims, numValues, isComplex, dataElement.get() );

    numRows    = dims[0];
    numColumns = dims[1];

    if ( matClass != MATIOStream::MAT_SPARSE_CLASS )
    {
        if ( MATIOStream::class2ScalarType( matClass ) == common::ScalarType::UNKNOWN )
        {
            COMMON_THROWEXCEPTION( "File " << mFile.getFileName() << " contains unsupported matrix class = " << matClass )
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

    if ( mScalarTypeData == common::ScalarType::PATTERN )
    {
        values.setSameValue( nnz, ValueType( 1 ) );   // set values with default value
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

    // build CSR storage that will be transposed

    CSRStorage<ValueType> csrStorage( dims[1], dims[0], std::move( ja ), std::move( ia ), std::move( values ) );

    csrStorage.assignTranspose( csrStorage );

    // transpose might change order, so we set back the diagonal elements as first entries

    csrStorage.setDiagonalFirst();

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
    else if ( MATIOStream::class2ScalarType( matClass ) == common::ScalarType::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File contains unsupported matrix class = " << matClass )
    }
    else
    {
        SCAI_LOG_INFO( logger, "Get dense<" << common::TypeTraits<ValueType>::stype << "> matrix "
                       << dims[0] << " x " << dims[1] )

        HArray<ValueType> values;

        offset += getArrayData( values, dataElementPtr + offset, nBytes - offset );

        if ( isComplex )
        {
            HArray<ValueType> imagValues;
            offset += getArrayData( imagValues, dataElementPtr + offset, nBytes - offset );
            buildComplex( values, imagValues );
        }

        // MATLAB stores it columnwise, so we transpose the data

        DenseStorage<ValueType> denseStorage( dims[1], dims[0], std::move( values ) );
        storage.assignTranspose( denseStorage );
    }

    SCAI_ASSERT_EQ_ERROR( offset, nBytes, "mismatch read bytes and size bytes, maybe COMPLEX" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readStorageImpl( MatrixStorage<ValueType>& storage )
{
    std::unique_ptr<char[]> dataElement;

    uint32_t nBytes = mFile.readDataElement( dataElement );

    getStorage( storage, dataElement.get(), nBytes );
    SCAI_LOG_INFO( logger, "readStorage: " << storage )
}

/* ------------------------------------------------------------------------------------ */
/*   Read grid vector data from file                                              */
/* ------------------------------------------------------------------------------------ */

void MatlabIO::readGridArray( _HArray& data, common::Grid& grid )
{
    // use the IO wrapper to call a typed version of readGrid

    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::readGridImpl( *this, data, grid );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readGridImpl( HArray<ValueType>& data, common::Grid& grid )
{
    std::unique_ptr<char[]> dataElement;

    uint32_t nBytes = mFile.readDataElement( dataElement );

    const char* elementPtr = dataElement.get();

    IndexType nDims;
    IndexType dims[SCAI_GRID_MAX_DIMENSION];
    IndexType nnz;
    bool      isComplex;

    MATIOStream::MATClass matClass;

    uint32_t offset = MATIOStream::getMatrixInfo( matClass, dims, SCAI_GRID_MAX_DIMENSION, nDims, nnz, isComplex, dataElement.get() );

    if ( matClass == MATIOStream::MAT_SPARSE_CLASS )
    {
        COMMON_THROWEXCEPTION( "File " << mFile.getFileName() << " contains sparse matrix, but not grid array" )
    }

    if ( MATIOStream::class2ScalarType( matClass ) == common::ScalarType::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "File " << mFile.getFileName() << " contains unsupported matrix class = " << matClass )
    }

    grid = common::Grid( nDims, dims );

    SCAI_ASSERT_LE_ERROR( offset, nBytes, "data element insufficient to read matrix info" )

    // now read the data

    offset += getArrayData( data, elementPtr + offset, nBytes - offset );

    SCAI_ASSERT_EQ_ERROR( data.size(), grid.size(), "serious mismatch" )

    if ( isComplex )
    {
        HArray<ValueType> imagValues;

        offset += getArrayData( imagValues, elementPtr + offset, nBytes - offset );

        buildComplex( data, imagValues );
    }

    HArray<ValueType> tmpData( grid.size() );
    changeMajor( tmpData, data, grid, false );
    data.swap( tmpData );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeStorage( const _MatrixStorage& storage )
{
    IOWrapper<MatlabIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readStorage( _MatrixStorage& storage )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatlabIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeArray( const hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( *this, array );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( *this, n, indexes, values );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readArray( hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( *this, array );
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatlabIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( *this, size, indexes, values );
}

/* --------------------------------------------------------------------------------- */

std::string MatlabIO::getMatrixFileSuffix() const
{
    return MatlabIO::createValue();
}

/* --------------------------------------------------------------------------------- */

std::string MatlabIO::getVectorFileSuffix() const
{
    return MatlabIO::createValue();
}

/* --------------------------------------------------------------------------------- */

#define SCAI_MATLAB_METHOD_INSTANTIATIONS( _type )           \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::writeArrayImpl(                           \
        const hmemo::HArray<_type>& array );                 \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::readArrayImpl(                            \
        hmemo::HArray<_type>& array );                       \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::writeSparseImpl(                          \
        const IndexType size,                                \
        const HArray<IndexType>& index,                      \
        const HArray<_type>& values );                       \
                                                             \
    template COMMON_DLL_IMPORTEXPORT                         \
    void MatlabIO::readSparseImpl(                           \
        IndexType& size,                                     \
        HArray<IndexType>& indexes,                          \
        HArray<_type>& values );                             \

SCAI_COMMON_LOOP( SCAI_MATLAB_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_MATLAB_METHOD_INSTANTIATIONS

#define SCAI_MATLAB_METHOD_INSTANTIATIONS( _type )      \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void MatlabIO::writeStorageImpl(                    \
        const MatrixStorage<_type>& storage );          \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void MatlabIO::readStorageImpl(                     \
        MatrixStorage<_type>& storage );                \

SCAI_COMMON_LOOP( SCAI_MATLAB_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_MATLAB_METHOD_INSTANTIATIONS

}  // lama

}  // scai


