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

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/io/IOStream.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Math.hpp>

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
    const common::scalar::ScalarType dataType )
{
    outFile << "%%MatrixMarket ";

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

        case common::scalar::INT:
        case common::scalar::LONG:
            outFile << "integer ";
            break;

        case common::scalar::PATTERN:
            outFile << "pattern ";
            break;

        default:
            COMMON_THROWEXCEPTION( *this << ": writeMMHeader: " << dataType << " unspported type" )
    }

    outFile << symmetry2str( symmetry ) << std::endl;

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

void MatrixMarketIO::readMMHeader(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    common::scalar::ScalarType& dataType,
    Symmetry& symmetry,
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

    if ( buffer == "real" )
    {
        dataType = common::scalar::FLOAT;
    }
    else if ( buffer == "double" ) 
    {
        dataType = common::scalar::DOUBLE;
    }
    else if ( buffer == "integer" ) 
    {
        dataType = common::scalar::INDEX_TYPE;
    }
    else if ( buffer == "complex" ) 
    {
        dataType = common::scalar::COMPLEX;
    }
    else if ( buffer == "pattern" ) 
    {
        dataType = common::scalar::PATTERN;
    }
    else 
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "Reading Matrix Market fille " << inFile.getFileName()
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

    common::scalar::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::scalar::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    bool      isVector   = true;
    IndexType numRows    = array.size();
    IndexType numColumns = 1;
    IndexType numValues  = -1;
    Symmetry  symmetry   = GENERAL;

    writeMMHeader( outFile, isVector, numRows, numColumns, numValues, symmetry, dataType );

    // output code runs only for host context

    ContextPtr host = Context::getHostPtr();

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( array, precData );

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixMarketIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName ) 
{
    Symmetry symmetry;
    common::scalar::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    ValueType val;
    std::string line;
    IOStream inFile( fileName, std::ios::in );
    readMMHeader( numRows, numColumns, numValues, mmType, symmetry, inFile );

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

        if ( mmType == common::scalar::PATTERN )
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

struct indexLess {

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
void MatrixMarketIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    SCAI_ASSERT_ERROR( mFileMode != BINARY, *this << ": Matrix market format can not be written binary" );

    COOStorage<ValueType> coo( storage );

    int numRows = coo.getNumRows();
    int numCols = coo.getNumColumns();

    // define empty array that will be swapped with the COOStorage

    LArray<IndexType> cooIA;
    LArray<IndexType> cooJA;
    LArray<ValueType> cooValues;

    coo.swap( cooIA, cooJA, cooValues );

    Symmetry symFlag = checkSymmetry( cooIA, cooJA, cooValues );

    SCAI_LOG_INFO( logger, "symmetry = " << symmetry2str( symFlag ) ) 

    if ( symFlag == SYMMETRIC || symFlag == HERMITIAN )
    {
        removeUpperTriangular( cooIA, cooJA, cooValues );
  
        SCAI_LOG_INFO( logger, "#values = " << cooIA.size() << ", due to symmetry " << symmetry2str( symFlag ) )
    }

    // Attention: indexing in MatrixMarket starts with 1 and not with 0 as in LAMA

    cooIA += 1;
    cooJA += 1;

    int numValues = cooIA.size();

    common::scalar::ScalarType dataType = mScalarTypeData;

    if ( dataType == common::scalar::INTERNAL )
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

    if ( dataType == common::scalar::PATTERN )
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

template<typename ValueType>
void MatrixMarketIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    Symmetry symmetry;
    common::scalar::ScalarType mmType;

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesFile;

    IOStream inFile( fileName, std::ios::in );

    readMMHeader( numRows, numColumns, numValuesFile, mmType, symmetry, inFile );

    SCAI_LOG_DEBUG( logger, "from header: nrows = " << numRows << ", ncols = " << numColumns << ", nnz = " << numValuesFile
                            << ", mmType = " << mmType << ", symmetry = " << symmetry )

    // check for consistency

    if ( common::isComplex( mmType ) && !common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        SCAI_LOG_WARN( logger, "Read matrix from Matrix Market file " << fileName 
                               << ": contains complex data but read in non-complex storage " << storage )
    }

    // use local arrays instead of heteregeneous arrays as we want ops on them

    LArray<IndexType> ia;   // row indexes, as ia in COO format
    LArray<IndexType> ja;   // col indexes, as ja in COO format
    LArray<ValueType> val;  // values as in COO format

    SCAI_LOG_DEBUG( logger, "read in" )
    
    if ( common::scalar::PATTERN == mmType )
    {
        inFile.readFormatted( ia, ja, numValuesFile );
        val.resize( numValuesFile );
        val = ValueType( 1 );
    }
    else
    {
        inFile.readFormatted( ia, ja, val, numValuesFile );
    }

    SCAI_LOG_DEBUG( logger, "read ia  : " << ia  )
    SCAI_LOG_DEBUG( logger, "read ja  : " << ja  )
    SCAI_LOG_DEBUG( logger, "read val : " << val )

    // MatrixMarket starts indexes always with one, so shift all row/col indexes

    ia -= 1;
    ja -= 1;

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

    // we shape the matrix by maximal appearing indexes

    int nrows = ia.max() + 1;
    int ncols = ja.max() + 1;

    SCAI_LOG_INFO( logger, "size from header: " << numRows << " x " << numColumns
                           << ", size by indexes: " << nrows << " x " << ncols )

    // specified size might be greater, but less is ERROR

    SCAI_ASSERT_GE( numRows, nrows, "found bigger row indexes than " << numRows )
    SCAI_ASSERT_GE( numColumns, ncols, "found bigger col indexes than " << numColumns )

    // take the COO arrays and build a COO storage that takes ownership of the data

    COOStorage<ValueType> coo( numRows, numColumns );

    coo.swap( ia, ja, val );

    storage = coo;
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
