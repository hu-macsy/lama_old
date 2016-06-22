/**
 * @file MatrixMarketIO.cpp
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
 * @brief Implementation of methods
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include "MatrixMarketIO.hpp"
#include "IOStream.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

#include <sstream>
#include <iomanip>

using namespace std;

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

static std::string MM_SUFFIX = ".mm";

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

void MatrixMarketIO::writeAt( std::ostream& stream ) const
{
    stream << "MatrixMarketIO ( ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

static void writeMMHeader(
    const bool& vector,
    const IndexType& numRows,
    const IndexType& numColumns,
    const IndexType& numValues,
    IOStream& outFile,
    const common::scalar::ScalarType& dataType )
{
    outFile << "%%matrixmarket ";

    if ( vector )
    {
        outFile << "vector array ";
    }
    else
    {
        outFile << "matrix coordinate ";
    }

    switch ( dataType )
    {
        case common::scalar::DOUBLE:
        case common::scalar::FLOAT:
        case common::scalar::LONG_DOUBLE:
            outFile << "real ";
            break;

        case common::scalar::COMPLEX:
        case common::scalar::DOUBLE_COMPLEX:
        case common::scalar::LONG_DOUBLE_COMPLEX:
            outFile << "complex ";
            break;

        case common::scalar::INDEX_TYPE:
            outFile << "integer ";
            break;

        case common::scalar::PATTERN:
            outFile << "pattern ";
            break;

        default:
            COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: " "unknown datatype." << dataType )
    }

    // TODO: Add support for symmetric matrices
    // currently we can only write non-symmetric
    outFile << "general" << std::endl;

    if ( vector )
    {
        outFile << numRows << " " << numColumns << std::endl;
    }
    else
    {
        outFile << numRows << " " << numColumns << " " << numValues << std::endl;
    }
}

/* --------------------------------------------------------------------------------- */

static void readMMHeader(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    bool& isPattern,
    bool& isSymmetric,
    IOStream& inFile )
{
    std::string buffer;
    // read %%MatrixMarket
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "%%matrixmarket" )
    {
        COMMON_THROWEXCEPTION( "Given file is no valid matrix market file, expected file to begin with %%MatrixMarket" )
    }

    // read object type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    // check if object type is valid in general
    if ( buffer != "matrix" && buffer != "vector" )
    {
        COMMON_THROWEXCEPTION( "Object type in the given matrix market file is invalid, should be matrix or vector" )
    }

    bool isVector = false;

    if ( buffer == "vector" )
    {
        isVector = true;
    }

    // read file type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    // checkif file type is valid in general
    if ( buffer != "coordinate" && buffer != "array" )
    {
        COMMON_THROWEXCEPTION( "Format type in the given matrix market file is invalid, should be coordinate or array" )
    }

    // read data type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "real" && buffer != "integer" && buffer != "complex" && buffer != "pattern" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be real, integer, complex or pattern" )
    }

    // TODO: allow to return other value types as well => check if the valid type is used
    if ( buffer == "pattern" )
    {
        isPattern = true;
    }
    else
    {
        isPattern = false;
    }

    // read symmetry
    std::getline( inFile, buffer, '\n' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "general" && buffer != "symmetric" && buffer != "skew-symmetric" && buffer != "hermitian" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be general, symmetric, skew-symmetric or hermitian" )
    }

    if ( buffer == "general" )
    {
        isSymmetric = false;
    }
    else
    {
        if ( buffer == "symmetric" )
        {
            isSymmetric = true;
        }
        else
        {
            // TODO: add support!
            COMMON_THROWEXCEPTION( "Symmetry options 'skew-symmetric' and 'hermitian' are currently not supported!" )
        }
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
    // TODO: vector correct here? should it be dense vs sparse?
    if ( !isVector )
    {
        bufferSS >> numValues;
    }
    else
    {
        numValues = numRows * numColumns;
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
static void addSymmetricEntries(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>& vals )
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
            wVals[offset ] = wVals[i];

            ++offset;
        }
    }

    ia.resize( offset );
    ja.resize( offset );
    vals.resize( offset );

    std::cout << "addSymmetricEntries: new size = " << offset << ", was " << numValues << std::endl;
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatrixMarketIO::logger, "FileIO.MatrixMarketIO" )

/* --------------------------------------------------------------------------------- */

bool MatrixMarketIO::isSupported( const bool binary ) const
{
    if ( binary )
    {
        return false; // binary is not supported
    }
    else
    {
        return true;  // formatted supported
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT_ERROR( !mBinary, "Matrix market format can not be written binary" );

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    common::scalar::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::scalar::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    IndexType numRows    = array.size();
    IndexType numColumns = 1;

    writeMMHeader( true, numRows, numColumns, -1, outFile, dataType );

    // output code runs only for host context

    ContextPtr host = Context::getHostPtr();

    ReadAccess<ValueType> dataRead( array, host );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        outFile << dataRead[i] << std::endl;
    }

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName ) 
{
    bool isSymmetric;
    bool isPattern;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    ValueType val;
    std::string line;
    IOStream inFile( fileName, std::ios::in );
    readMMHeader( numRows, numColumns, numValues, isPattern, isSymmetric, inFile );

    if ( numColumns != 1 )
    {
        SCAI_LOG_WARN( logger, "reading vector from mtx file, #columns = " << numColumns << ", ignored" )
    }
 
    // Note: we ignore number of columns here and make a vector of size numRows x numColumns

    WriteOnlyAccess<ValueType> vector( array, numValues );

    IndexType i;
    ValueType* vPtr = vector.get();

    for ( int l = 0; l < numValues && !inFile.eof(); ++l )
    {
        std::getline( inFile, line );
        std::istringstream reader( line );

        if ( isPattern )
        {
            reader >> i;
            val = 1.0;
            i--;
        }
        else
        {
            reader >> val;
            i = l;
        }

        vPtr[i] = val;
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

template<typename ValueType>
void MatrixMarketIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    SCAI_ASSERT( !mBinary, "Binary mode not supported for MatrixMarketIO" )

    COOStorage<ValueType> coo( storage );

    int numRows = coo.getNumRows();
    int numCols = coo.getNumColumns();

    // define empty array that will be swapped with the COOStorage

    LArray<IndexType> cooIA;
    LArray<IndexType> cooJA;
    LArray<ValueType> cooValues;

    coo.swap( cooIA, cooJA, cooValues );

    // Attention: indexing in MatrixMarket starts with 1 and not with 0 as in LAMA

    cooIA += 1;
    cooJA += 1;

    int numValues = cooIA.size();

    common::scalar::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::scalar::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    writeMMHeader( false, numRows, numCols, numValues, outFile, dataType );

    // output code runs only for host context

    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( cooIA, host );
    ReadAccess<IndexType> ja( cooJA, host );
    ReadAccess<ValueType> data( cooValues, host );

    int precIndex = 0;
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex, cooValues, precData );

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    bool isSymmetric;
    bool isPattern;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesFile;

    IOStream inFile( fileName, std::ios::in );

    readMMHeader( numRows, numColumns, numValuesFile, isPattern, isSymmetric, inFile );

    SCAI_LOG_DEBUG( logger, "from header: nrows = " << numRows << ", ncols = " << numColumns << ", nnz = " << numValuesFile
                            << ", isPattern = " << isPattern << ", isSymmetric = " << isSymmetric )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    LArray<IndexType> ia;   // row indexes, as ia in COO format
    LArray<IndexType> ja;   // col indexes, as ja in COO format
    LArray<ValueType> val;  // values as in COO format

    SCAI_LOG_DEBUG( logger, "read in" )
    
    inFile.readFormatted( ia, ja, val, numValuesFile );

    SCAI_LOG_DEBUG( logger, "read ia  : " << ia  )
    SCAI_LOG_DEBUG( logger, "read ja  : " << ja  )
    SCAI_LOG_DEBUG( logger, "read val : " << val )

    // MatrixMarket starts indexes always with one, so shift all row/col indexes

    ia -= 1;
    ja -= 1;

    // double symmetric entries

    if ( isSymmetric )
    {
        // add symmetric entries, no check for doubles

        addSymmetricEntries( ia, ja, val );

        SCAI_LOG_DEBUG( logger, "sym ia  : " << ia  )
        SCAI_LOG_DEBUG( logger, "sym ja  : " << ja  )
        SCAI_LOG_DEBUG( logger, "sym val : " << val )
    }

    // check if there is more data in the file tht should not be there

    std::string line;

    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': invalid file, unread lines" )
    }

    inFile.close();

    // we shape the matrix by maximal appearing indexes

    int nrows = ia.max() + 1;
    int ncols = ja.max() + 1;

    SCAI_LOG_INFO( logger, "size from header: " << numRows << " x " << numColumns
                           << ", size by indexes: " << nrows << " x " << ncols )

    // specified size might be greater, but less is ERROR

    SCAI_ASSERT_GE( numRows, nrows, "found bigger row indexes than " << numRows )
    SCAI_ASSERT_GE( numColumns, ncols, "found bigger col indexes than " << numColumns )

    COOStorage<ValueType> coo( numRows, numColumns, ia, ja, val );

    storage = coo;
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
