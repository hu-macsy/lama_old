/**
 * @file MatrixMarketIO.cpp
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
 * @brief Implementation of methods to read/write Matrix Market files
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include <scai/lama/io/MatrixMarketIO.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Grid.hpp>

#include <scai/tracing.hpp>

#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

#define MM_SUFFIX ".mtx"

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

std::string MatrixMarketIO::getVectorFileSuffix() const
{
    return MM_SUFFIX;
}

std::string MatrixMarketIO::getMatrixFileSuffix() const
{
    return MM_SUFFIX;
}

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* MatrixMarketIO::create()
{
    return new MatrixMarketIO();
}

std::string MatrixMarketIO::createValue()
{
    return MM_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool MatrixMarketIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeAt( std::ostream& stream ) const
{
    stream << "MatrixMarketIO ( ";
    stream << "suffix = " << MM_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

const char* MatrixMarketIO::symmetry2str( const Symmetry symmetry )
{
    switch ( symmetry )
    {
        case GENERAL:
            return "general";

        case SYMMETRIC:
            return "symmetric";

        case HERMITIAN:
            return "hermitian";

        case SKEW_SYMMETRIC:
            return "skew-symmetric";

        default:
            return "unknown";
    }
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeMMHeader(
    IOStream& outFile,
    const bool vector,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const Symmetry symmetry,
    const common::ScalarType dataType )
{
    outFile << "%%MatrixMarket ";

    if ( vector )
    {
        outFile << "vector ";
    }
    else
    {
        outFile << "matrix ";
    }

    if ( numValues == invalidIndex )
    {
        outFile << "array ";
    }
    else
    {
        outFile << "coordinate ";
    }

    switch ( dataType )
    {
        case common::ScalarType::DOUBLE:
        case common::ScalarType::FLOAT:
        case common::ScalarType::LONG_DOUBLE:
            outFile << "real ";
            break;

        case common::ScalarType::COMPLEX:
        case common::ScalarType::DOUBLE_COMPLEX:
        case common::ScalarType::LONG_DOUBLE_COMPLEX:
            outFile << "complex ";
            break;

        case common::ScalarType::INT:
        case common::ScalarType::LONG:
            outFile << "integer ";
            break;

        case common::ScalarType::PATTERN:
            outFile << "pattern ";
            break;

        default:
            COMMON_THROWEXCEPTION( *this << ": writeMMHeader: " << dataType << " unspported type" )
    }

    outFile << symmetry2str( symmetry ) << std::endl;

    // array ( numValues == invalidIndex ) has no entry numValues, only coordinate

    if ( numValues == invalidIndex )
    {
        outFile << numRows << " " << numColumns << std::endl;
    }
    else
    {
        outFile << numRows << " " << numColumns << " " << numValues << std::endl;
    }
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readMMHeader(
    IOStream& inFile,
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    common::ScalarType& dataType,
    bool& isVector,
    Symmetry& symmetry )
{
    const std::string& fileName = inFile.getFileName();

    std::string buffer;

    // Header: %%MatrixMarket ...

    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "%%matrixmarket" )
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "File " << fileName << " is no valid matrix market file"
                             << ", expected file to begin with %%MatrixMarket" )
    }

    // read object type matrix | vector

    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    // check if object type is valid in general

    if ( buffer == "matrix" )
    {
        isVector = false;
    }
    else if ( buffer == "vector" )
    {
        isVector = true;
    }
    else
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "Reading Matrix Market file " << fileName << ": object type " << buffer
                             << " illegal, must be matrix or vector" )
    }

    // read file type coordinate | array

    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    bool isCoordinate = false;

    // check if file type is valid in general

    if ( buffer == "coordinate" )
    {
        isCoordinate = true;
    }
    else if ( buffer == "array" )
    {
        isCoordinate = false;
    }
    else
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "Reading Matrix Market file " << fileName << ": format type " << buffer
                             << " illegal, must be array or coordinate" )
    }

    // read data type

    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer == "real" )
    {
        dataType = common::ScalarType::FLOAT;
    }
    else if ( buffer == "double" )
    {
        dataType = common::ScalarType::DOUBLE;
    }
    else if ( buffer == "integer" )
    {
        dataType = common::ScalarType::INDEX_TYPE;
    }
    else if ( buffer == "complex" )
    {
        dataType = common::ScalarType::COMPLEX;
    }
    else if ( buffer == "pattern" )
    {
        dataType = common::ScalarType::PATTERN;
    }
    else
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "Reading Matrix Market fille " << fileName
                             << ": data type field = " << buffer << " is illegal"
                             << ", should be real, double, integer, complex, pattern" )
    }

    // read symmetry

    std::getline( inFile, buffer, '\n' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer == "general" )
    {
        symmetry = GENERAL;
    }
    else if ( buffer == "symmetric" )
    {
        symmetry = SYMMETRIC;
    }
    else if ( buffer == "hermitian" )
    {
        symmetry = HERMITIAN;
    }
    else if ( buffer == "skew-symmetric" )
    {
        symmetry = SKEW_SYMMETRIC;
    }
    else
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "Reading Matrix Market fille " << inFile.getFileName()
                             << ": symmetry = " << buffer << " is illegal"
                             << ", should be general, symmetric, skew-symmetric or hermitian" )
    }

    // skip further comment lines

    bool skip = true;

    do
    {
        std::getline( inFile, buffer, '\n' );

        if ( inFile.fail() )
        {
            // no further line, that is serious

            COMMON_THROWEXCEPTION( "line with matrix / vector sizes not found" )
        }

        if ( buffer.size() > 0 )
        {
            skip = buffer[0] == '%';
        }
    }
    while ( skip );

    std::stringstream bufferSS( buffer );

    bufferSS >> numRows;
    bufferSS >> numColumns;

    if ( isCoordinate )
    {
        bufferSS >> numValues;
    }
    else
    {
        numValues = invalidIndex;   // stands for array
    }

    if ( isVector && numColumns != 1 )
    {
        SCAI_LOG_WARN( logger, "MatrixMarket file " << fileName << ": vector, but ncol = " << numColumns )
    }

    SCAI_LOG_INFO( logger, "read MM header: size = " << numRows << " x " << numColumns
                   << ", #nnz = " << numValues << ", type = " << dataType
                   << ", symmetry = " << symmetry )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::addSymmetricEntries(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>& vals,
    bool conjFlag )
{
    IndexType numValues = ia.size();

    SCAI_ASSERT_EQUAL( numValues, ja.size(), "size mismatch" );
    SCAI_ASSERT_EQUAL( numValues, vals.size(), "size mismatch" );

    // make sure that sufficient memory is available

    ContextPtr host = Context::getHostPtr();

    // reserve of memory gurantees that we can already use it

    ia.reserve( host, 2 * numValues );
    ja.reserve( host, 2 * numValues );
    vals.reserve( host, 2 * numValues );

    WriteAccess<IndexType> wIA( ia, host );
    WriteAccess<IndexType> wJA( ja, host );
    WriteAccess<ValueType> wVals( vals, host );

    IndexType offset = numValues;   // index for added entries

    for ( IndexType i = 0; i < numValues; ++i )
    {
        if ( wIA[i] != wJA[i] )
        {
            wIA[ offset ] = wJA[i];
            wJA[ offset ] = wIA[i];

            if ( conjFlag )
            {
                wVals[offset ] = common::Math::conj( wVals[i] );
            }
            else
            {
                wVals[offset ] = wVals[i];
            }

            ++offset;
        }
    }

    ia.resize( offset );
    ja.resize( offset );
    vals.resize( offset );

    SCAI_LOG_INFO( logger, "addSymmetricEntries: new size = " << offset << ", was " << numValues )
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatrixMarketIO::logger, "FileIO.MatrixMarketIO" )

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_LOG_INFO( logger, *this << ": write array " << array << " to " << fileName );

    SCAI_ASSERT_ERROR( mFileMode != BINARY, *this << ": Matrix market format can not be written binary" );

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    common::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::ScalarType::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    bool      isVector   = true;
    IndexType numRows    = array.size();
    IndexType numColumns = 1;
    IndexType numValues  = invalidIndex;
    Symmetry  symmetry   = GENERAL;

    writeMMHeader( outFile, isVector, numRows, numColumns, numValues, symmetry, dataType );

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( array, precData );

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::writeSparseImpl(
    const IndexType size,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values,
    const std::string& fileName )
{
    SCAI_ASSERT_EQ_ERROR( indexes.size(), values.size(), "size mismatch for indexes/values in sparse array" )

    IndexType numValues  = indexes.size();

    SCAI_LOG_INFO( logger, *this << ": write sparse to " << fileName );

    SCAI_ASSERT_ERROR( mFileMode != BINARY, *this << ": Matrix market format can not be written binary" );

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    common::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::ScalarType::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    bool      isVector   = true;
    IndexType numRows    = size;
    IndexType numColumns = 1;
    Symmetry  symmetry   = GENERAL;    // does not matter for Vector

    writeMMHeader( outFile, isVector, numRows, numColumns, numValues, symmetry, dataType );

    int precIndexes = getDataPrecision( indexes.getValueType() );
    int precValues  = getDataPrecision( values.getValueType() );

    HArray<IndexType> indexes1;

    HArrayUtils::compute<IndexType>( indexes1, indexes, common::BinaryOp::ADD, 1 );

    outFile.writeFormatted( indexes1, precIndexes, values, precValues );

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readArrayInfo( IndexType& size, const std::string& fileName )
{
    SCAI_REGION( "IO.MM.readArrayInfo" )

    Symmetry symmetry;
    common::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;
    bool      isVector;

    IOStream inFile( fileName, std::ios::in );

    readMMHeader( inFile, numRows, numColumns, numValues, mmType, isVector, symmetry );

    if ( !isVector )
    {
        SCAI_LOG_WARN( logger, "Matrix Market file " << fileName << ", contains matrix and not vector" )
    }

    size = numRows * numColumns;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readVectorCoordinates( 
    HArray<IndexType>& indexes, 
    HArray<ValueType>& values,
    IOStream& inFile, 
    const IndexType numValues, 
    const bool isVector,
    common::ScalarType mmType )
{
    // todo: pattern for vector

    SCAI_ASSERT_ERROR( mmType != common::ScalarType::PATTERN, "pattern not handled yet" )

    if ( isVector )
    {
        // only one position in each line

        inFile.readFormatted( indexes, values, numValues );
    }
    else
    {
        // row, col position in each line

        HArray<IndexType> colDummy;
        inFile.readFormatted( indexes, colDummy, values, numValues );
    }

    HArrayUtils::setScalar<IndexType>( indexes, 1, common::BinaryOp::SUB );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName,
    const IndexType first,
    const IndexType n )
{
    SCAI_REGION( "IO.MM.readArray" )

    Symmetry symmetry;
    common::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;
    bool      isVector;

    std::string line;
    IOStream inFile( fileName, std::ios::in );
    readMMHeader( inFile, numRows, numColumns, numValues, mmType, isVector, symmetry );

    IndexType size = numRows * numColumns;

    if ( numColumns != 1 )
    {
        SCAI_LOG_WARN( logger, "reading vector from mtx file, #columns = " << numColumns << ", ignored" )
    }

    // Note: we ignore number of columns here and make a vector of size numRows x numColumns

    if ( ! common::Utils::validIndex( first, size ) )
    {
        array.clear();
        return;
    }

    IndexType nEntries = n;

    if ( n == invalidIndex )
    {
        nEntries = size - first;
    }
    else
    {
        SCAI_ASSERT_LE_ERROR( first + n, size, "array block size " << n << " invalid" )
    }

    if ( numValues != invalidIndex )
    {
        // so we have coordinate format

        HArray<IndexType> indexes;
        HArray<ValueType> values;

        readVectorCoordinates( indexes, values, inFile, numValues, isVector, mmType );
 
        HArrayUtils::buildDenseArray( array, size, values, indexes, ValueType( 0 ) );
    }
    else
    {
        // so we have dense array format

        inFile.readFormatted( array, size );
    }

    if ( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': invalid file, contains to many elements." )
    }

    inFile.close();
    SCAI_LOG_INFO( logger, "read array " << numRows )

    if ( nEntries != size )
    {
        hmemo::HArray<ValueType> block( nEntries );
        hmemo::ContextPtr ctx = hmemo::Context::getHostPtr();
        SCAI_LOG_DEBUG( logger, "read block first = " << first << ", n = " << nEntries << " from array " << array )

        IndexType inc = 1;
        utilskernel::HArrayUtils::setArraySection( block, 0, inc, array, first, inc, nEntries, common::BinaryOp::COPY, ctx );

        array.swap( block );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const std::string& fileName )
{
    SCAI_REGION( "IO.MM.readSparse" )

    Symmetry symmetry;
    common::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;
    bool      isVector;

    std::string line;
    IOStream inFile( fileName, std::ios::in );
    readMMHeader( inFile, numRows, numColumns, numValues, mmType, isVector, symmetry );

    size = numRows * numColumns;

    if ( numColumns != 1 )
    {
        SCAI_LOG_WARN( logger, "reading vector from mtx file, #columns = " << numColumns << ", ignored" )
    }

    if ( numValues != invalidIndex )
    {
        // so we have coordinate format

        readVectorCoordinates( indexes, values, inFile, numValues, isVector, mmType );
    }
    else
    {
        // so we have dense array format

        HArray<ValueType> denseArray;
        inFile.readFormatted( denseArray, size );
        ValueType zero = 0;
        HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
    }

    if ( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': invalid file, contains to many elements." )
    }

    inFile.close();

    SCAI_LOG_INFO( logger, "read array " << numRows )
}

/* --------------------------------------------------------------------------------- */

struct indexLess
{

    const IndexType* ia;
    const IndexType* ja;

    bool operator()( int pos1, int pos2 )
    {
        return    ( ia[pos1] < ia[pos2] )
                  || ( ia[pos1] == ia[pos2] && ja[pos1] < ja[pos2] );
    }
};

static void sortIJ( IndexType perm[], const IndexType ia[], const IndexType ja[], IndexType N )
{
    for ( IndexType i = 0; i < N; ++ i )
    {
        perm[i] = i;
    }

    indexLess cmp;

    cmp.ia = ia;
    cmp.ja = ja;

    // sort using a custom function object

    std::sort( perm, perm + N, cmp );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
MatrixMarketIO::Symmetry MatrixMarketIO::checkSymmetry( const HArray<IndexType>& cooIA, const HArray<IndexType>& cooJA, const HArray<ValueType>& cooValues )
{
    SCAI_REGION( "IO.MM.checkSymmetry" )

    IndexType n = cooIA.size();

    bool isSym  = true;
    bool isHerm = true;

    HArray<IndexType> rank1;   // sorted IA, JA
    HArray<IndexType> rank2;   // sorted JA, IA

    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( cooIA, host );
    ReadAccess<IndexType> ja( cooJA, host );
    ReadAccess<ValueType> values( cooValues, host );

    WriteOnlyAccess<IndexType> perm1( rank1, host, n );
    WriteOnlyAccess<IndexType> perm2( rank2, host, n );

    sortIJ( perm1, ia, ja, n );
    sortIJ( perm2, ja, ia, n );

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType i1 = ia[perm1[i]];
        IndexType i2 = ja[perm2[i]];
        IndexType j1 = ja[perm1[i]];
        IndexType j2 = ia[perm2[i]];

        SCAI_LOG_TRACE( logger, "Check: pos1 = " << perm1[i] << ": ( " << i1 << ", " << j1
                        << " ), pos2 = " << perm2[i] << ": ( " << i2 << ", " << j2 << " )" )

        if ( i1 != i2 || j1 != j2 )
        {
            isSym  = false;
            isHerm = false;
            SCAI_LOG_DEBUG( logger, "entry (" << i1 << ", " << j1 << ") available, "
                            << "but not entry (" << j1 << ", " << i1 << ")" )
            break;
        }

        // we have found entry( i1, j1 ) and entry( j1, i1 )

        if ( i1 <= j1 )
        {
            // further check only for lower triangular part

            continue;
        }

        ValueType v1 = values[perm1[i]];
        ValueType v2 = values[perm2[i]];

        SCAI_LOG_TRACE( logger, "compare mirrored values " << v1 << " " << v2 )

        if ( v1 != v2 )
        {
            isSym = false;
        }

        if ( common::Math::conj( v1 ) != v2 )
        {
            isHerm = false;
        }

        if ( !isSym && !isHerm )
        {
            break;
        }
    }

    if ( isSym )
    {
        return SYMMETRIC;
    }
    else if ( isHerm )
    {
        return HERMITIAN;
    }
    else
    {
        return GENERAL;
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
static void removeUpperTriangular( HArray<IndexType>& cooIA, HArray<IndexType>& cooJA, HArray<ValueType>& cooValues )
{
    IndexType n = cooIA.size();

    IndexType k = 0;

    // take only entries( i, j ) with i >= j

    {
        ContextPtr host = Context::getHostPtr();

        WriteAccess<IndexType> ia( cooIA, host );
        WriteAccess<IndexType> ja( cooJA, host );
        WriteAccess<ValueType> values( cooValues, host );

        for ( IndexType pos = 0; pos < n; ++pos )
        {
            if ( ia[pos] >= ja[pos] )
            {
                ia[k] = ia[pos];
                ja[k] = ja[pos];
                values[k] = values[pos];
                k++;
            }
        }
    }

    // WriteAccess should have been free otherwise resize might not be allowed

    cooIA.resize( k );
    cooJA.resize( k );
    cooValues.resize( k );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::writeDenseMatrix(
    const DenseStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_REGION( "IO.MM.writeDense" )

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    common::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::ScalarType::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    bool      isVector   = false;
    IndexType numRows    = storage.getNumRows();
    IndexType numColumns = storage.getNumColumns();
    IndexType numValues  = invalidIndex;
    Symmetry  symmetry   = GENERAL;

    SCAI_LOG_INFO( logger, "write dense matrix " << numRows << " x " << numColumns )

    writeMMHeader( outFile, isVector, numRows, numColumns, numValues, symmetry, dataType );

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    // matrix market writes data column wise, so transpose it to get it in correct order

    DenseStorage<ValueType> transposedDenseStorage;
    transposedDenseStorage.assignTranspose( storage );

    outFile.writeFormatted( transposedDenseStorage.getData(), precData );

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    if ( storage.getFormat() == Format::DENSE )
    {
        const DenseStorage<ValueType>& denseStorage = reinterpret_cast<const DenseStorage<ValueType>&>( storage );
        writeDenseMatrix( denseStorage, fileName );
        return;
    }

    SCAI_REGION( "IO.MM.writeCOO" )

    SCAI_ASSERT_ERROR( mFileMode != BINARY, *this << ": Matrix market format can not be written binary" );

    auto coo = convert<COOStorage<ValueType>>( storage );   // converts (any) storage to COO format

    IndexType numRows = coo.getNumRows();
    IndexType numCols = coo.getNumColumns();

    // define arrays that will contain the COO daa

    HArray<IndexType> cooIA;       // = coo.getIA()
    HArray<IndexType> cooJA;       // = coo.getJA()
    HArray<ValueType> cooValues;   // = coo.getValues

    // use splitUp to avoid copies of the data 

    coo.splitUp( numRows, numCols, cooIA, cooJA, cooValues );

    Symmetry symFlag = checkSymmetry( cooIA, cooJA, cooValues );

    SCAI_LOG_INFO( logger, "symmetry = " << symmetry2str( symFlag ) )

    if ( symFlag == SYMMETRIC || symFlag == HERMITIAN )
    {
        removeUpperTriangular( cooIA, cooJA, cooValues );

        SCAI_LOG_INFO( logger, "#values = " << cooIA.size() << ", due to symmetry " << symmetry2str( symFlag ) )
    }

    // Attention: indexing in MatrixMarket starts with 1 and not with 0 as in LAMA

    HArrayUtils::setScalar<IndexType>( cooIA, 1, common::BinaryOp::ADD );
    HArrayUtils::setScalar<IndexType>( cooJA, 1, common::BinaryOp::ADD );

    int numValues = cooIA.size();

    common::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::ScalarType::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    // If file type is no more complex, the storage cannot be HERMITIAN any more

    if ( symFlag == HERMITIAN && !isComplex( dataType ) )
    {
        symFlag = SYMMETRIC;
    }

    bool isVector = false;

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    writeMMHeader( outFile, isVector, numRows, numCols, numValues, symFlag, dataType );

    // output code runs only for host context

    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( cooIA, host );
    ReadAccess<IndexType> ja( cooJA, host );
    ReadAccess<ValueType> data( cooValues, host );

    int precIndex = 0;

    if ( dataType == common::ScalarType::PATTERN )
    {
        outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex );
    }
    else
    {
        int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

        outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex, cooValues, precData );
    }

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readStorageInfo(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    const std::string& fileName )

{
    SCAI_REGION( "IO.MM.readStorageInfo" )

    Symmetry symmetry;
    common::ScalarType mmType;
    bool isVector;

    IndexType numValuesFile;

    IOStream inFile( fileName, std::ios::in );

    readMMHeader( inFile, numRows, numColumns, numValuesFile, mmType, isVector, symmetry );

    if ( numValuesFile == invalidIndex )
    {
        numValues = numRows * numColumns;
    }
    else
    {
        numValues = numValuesFile;

        if ( symmetry )
        {
            // in case of symmetry we assume one entry for each diagonal element and double the non-diagonal elements

            SCAI_ASSERT_EQUAL( numRows, numColumns, "symmetry only possible for square matrices" )

            numValues = numRows + 2 * ( numValuesFile - numRows );
        }
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readMMArray(
    IOStream& inFile,
    HArray<ValueType>& data,
    const IndexType numRows,
    const IndexType numColumns,
    const Symmetry symmetry )
{
    SCAI_REGION( "IO.MM.readDense" )

    const std::string& fileName = inFile.getFileName();

    if ( symmetry != GENERAL )
    {
        SCAI_ASSERT_EQ_ERROR( numRows, numColumns, "symmetric data only for square matrices" )
    }

    if ( symmetry == SKEW_SYMMETRIC )
    {
        COMMON_THROWEXCEPTION( "skew symmentric not supported yet" )
    }

    WriteOnlyAccess<ValueType> wData( data, numRows * numColumns );

    // values in input file are column-major order

    std::string line;   // used for reading lines of file

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        for ( IndexType i = 0; i < numRows; ++i )
        {
            // symmemtric data has no upper triangular, i < j

            if ( symmetry == SYMMETRIC && i < j )
            {
                wData[i * numColumns + j] = wData[ j * numRows + i];
                continue;
            }

            if ( symmetry == HERMITIAN && i < j )
            {
                wData[i * numColumns + j] = common::Math::conj( wData[ j * numRows + i] );
                continue;
            }

            if ( inFile.eof() )
            {
                COMMON_THROWEXCEPTION( "failed to read further entry" )
            }

            std::getline( inFile, line );
            std::istringstream reader( line );

            ValueType val;

            reader >> val;
            wData[i * numColumns + j] = val;
        }
    }

    if ( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': invalid file, contains to many elements." )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName,
    const IndexType firstRow,
    const IndexType nRows )
{
    SCAI_REGION( "IO.MM.readStorage" )

    Symmetry symmetry;
    common::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesFile;
    bool isVector;

    IOStream inFile( fileName, std::ios::in );

    readMMHeader( inFile, numRows, numColumns, numValuesFile, mmType, isVector, symmetry );

    SCAI_LOG_DEBUG( logger, "from header: nrows = " << numRows << ", ncols = " << numColumns << ", nnz = " << numValuesFile
                    << ", mmType = " << mmType << ", symmetry = " << symmetry )

    // check for consistency

    if ( common::isComplex( mmType ) && !common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        SCAI_LOG_WARN( logger, "Read matrix from Matrix Market file " << fileName
                       << ": contains complex data but read in non-complex storage " << storage )
    }

    if ( numValuesFile == invalidIndex )
    {
        HArray<ValueType> data;

        readMMArray( inFile, data, numRows, numColumns, symmetry );

        DenseStorage<ValueType> denseStorage( numRows, numColumns, std::move( data ) );

        if ( firstRow == 0 && nRows == invalidIndex )
        {
            storage = denseStorage;
        }
        else
        {
            denseStorage.copyBlockTo( storage, firstRow, nRows );
        }

        return;
    }

    // use local arrays instead of heteregeneous arrays as we want ops on them

    HArray<IndexType> ia;   // row indexes, as ia in COO format
    HArray<IndexType> ja;   // col indexes, as ja in COO format
    HArray<ValueType> val;  // values as in COO format

    SCAI_LOG_DEBUG( logger, "read in" )

    if ( common::ScalarType::PATTERN == mmType )
    {
        // to be consistent with other FileIO handlers throw an exception if not Pattern expected

        SCAI_ASSERT_EQ_ERROR( common::ScalarType::PATTERN, mScalarTypeData, "File " << fileName << " has only matrix pattern" )

        inFile.readFormatted( ia, ja, numValuesFile );
        val.setSameValue( numValuesFile, ValueType( 1 ) );
    }
    else
    {
        inFile.readFormatted( ia, ja, val, numValuesFile );
    }

    SCAI_LOG_DEBUG( logger, "read ia  : " << ia  )
    SCAI_LOG_DEBUG( logger, "read ja  : " << ja  )
    SCAI_LOG_DEBUG( logger, "read val : " << val )

    // MatrixMarket starts indexes always with one, so shift all row/col indexes

    HArrayUtils::setScalar<IndexType>( ia, 1, common::BinaryOp::SUB );
    HArrayUtils::setScalar<IndexType>( ja, 1, common::BinaryOp::SUB );

    // double symmetric entries

    if ( symmetry == SYMMETRIC )
    {
        // add symmetric entries, no check for doubles

        addSymmetricEntries( ia, ja, val, false );
    }
    else if ( symmetry == HERMITIAN )
    {
        // add hermitian entries, no check for doubles

        addSymmetricEntries( ia, ja, val, true );
    }
    else if ( symmetry == SKEW_SYMMETRIC )
    {
        // skew-symmetric not supported

        SCAI_THROWEXCEPTION( common::IOException, "Matrix Market file " << fileName
                             << ": skew-symmetric not supported yet" )
    }

    // close and check if there is more data in the file tht should not be there

    inFile.closeCheck();

    // we shape the matrix by maximal appearing indexes, apply max only if size > 0

    IndexType nrows = ia.size() ? HArrayUtils::max( ia ) + 1 : 0;
    IndexType ncols = ja.size() ? HArrayUtils::max( ja ) + 1 : 0;

    SCAI_LOG_INFO( logger, "size from header: " << numRows << " x " << numColumns
                   << ", size by indexes: " << nrows << " x " << ncols )

    // specified size might be greater, but less is ERROR

    SCAI_ASSERT_GE( numRows, nrows, "found bigger row indexes than " << numRows )
    SCAI_ASSERT_GE( numColumns, ncols, "found bigger col indexes than " << numColumns )

    // take the COO arrays and build a COO storage that takes ownership of the data

    COOStorage<ValueType> coo( numRows, numColumns, std::move( ia ), std::move( ja ), std::move( val ) );

    if ( firstRow == 0 && nRows == invalidIndex )
    {
        storage = coo;
    }
    else
    {
        coo.copyBlockTo( storage, firstRow, nRows );
    }
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data, outputFileName );
}

void MatrixMarketIO::readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    readArray( data, inputFileName );

    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "MatrixMarket does not support multidimensional array, take default shape " << grid )
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeStorage( const _MatrixStorage& storage, const std::string& fileName )
{
    IOWrapper<MatrixMarketIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( ( MatrixMarketIO& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readStorage(
    _MatrixStorage& storage,
    const std::string& fileName,
    const IndexType offsetRow,
    const IndexType nRows )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatrixMarketIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( ( MatrixMarketIO& ) *this, storage, fileName, offsetRow, nRows );
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeArray( const hmemo::_HArray& array, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( ( MatrixMarketIO& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( ( MatrixMarketIO& ) *this, n, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readArray( hmemo::_HArray& array, const std::string& fileName, const IndexType offset, const IndexType n )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( ( MatrixMarketIO& ) *this, array, fileName, offset, n );
}

/* --------------------------------------------------------------------------------- */

void MatrixMarketIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( ( MatrixMarketIO& ) *this, size, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
