/**
 * @file TextIO.cpp
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
 * @brief Implementation of methods for FileIO class TextIO
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include <scai/lama/io/TextIO.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/lama/storage/COOStorage.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Grid.hpp>
#include <scai/common/exception/IOException.hpp>

#include <sstream>

#define MATLAB_SUFFIX ".txt"

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

FileIO* TextIO::create()
{
    return new TextIO();
}

std::string TextIO::createValue()
{
    return MATLAB_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool TextIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeAt( std::ostream& stream ) const
{
    stream << "TextIO ( suffix = " << MATLAB_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/* --------------------------------------------------------------------------------- */

/** Method to count number of lines of a text file and the maximal number of entries in one line
 *
 *  @param[out]  nLines is number of lines the file has
 *  @param[out]  nEntries is maximal number of entries
 *  @param[in]   fileName is the name of the file
 *
 *  Note: it might be possible that one line contains less than 'nEntries' entries
 */
void TextIO::checkTextFile( IndexType& nLines, IndexType& nEntries, const char* fileName )
{
    nLines   = 0;
    nEntries = 0;

    std::ifstream infile( fileName, std::ios::in );

    if ( infile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    std::string line;
    std::vector<std::string> tokens;

    while ( std::getline( infile, line ) )
    {
        ++nLines;

        common::Settings::tokenize( tokens, line );

        IndexType nTokens = tokens.size();

        if ( nTokens > nEntries )
        {
            nEntries = nTokens;
            // LOG_DEBUG: cout << "max tokens = " << nEntries << " at line " << nLines << endl;
        }
    }

    SCAI_LOG_INFO( logger, "checkTextFile " << fileName << ": #lines = " << nLines << ", #entries = " << nEntries )
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( TextIO::logger, "FileIO.TextIO" )

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    IOStream outFile( fileName, std::ios::out );

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( array, precData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeSparseImpl(
    const IndexType size,
    const hmemo::HArray<IndexType>& indexes,
    const hmemo::HArray<ValueType>& values,
    const std::string& fileName )
{
    // Note: size is ignored for TextIO of sparse vector

    if ( true )
    {
        HArray<ValueType> denseArray;
        ValueType zero = 0;
        HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );
        writeArrayImpl( denseArray, fileName );
    }
    else
    {
        // Better solution to write sparse data, but not unique when reading, no header

        SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

        IOStream outFile( fileName, std::ios::out );
    
        int precIndexes  = getDataPrecision( indexes.getValueType() );
        int precData     = getDataPrecision( values.getValueType() );

        outFile.writeFormatted( indexes, precIndexes, values, precData );
    }
}

/* --------------------------------------------------------------------------------- */

void TextIO::readArrayInfo( IndexType& size, const std::string& fileName )
{
    IndexType nEntries;   // dummy variable needed for checkTextFile

    // each array entry in one line, so count the number of lines in the file

    checkTextFile( size, nEntries, fileName.c_str() );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName,
    const IndexType first,
    const IndexType n )
{
    IndexType size;   // number of lines, size of array
    IndexType k;      // number of entries in one line

    checkTextFile( size, k, fileName.c_str() );

    SCAI_LOG_INFO( logger, "File : " << fileName << ", #lines = " << size << ", #entries = " << k )

    SCAI_ASSERT_LE( k, 2, "#entries/row in file " << fileName << " must not excced 2" )

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
        // give useful error message as this is a typical error if wrong file is specified

        SCAI_ASSERT_LE_ERROR( first + n, size, "Read array block( offset = " << first << ", n = " << nEntries << ") failed: "
                              << "size of array in file " << fileName << " is " << size )
    }

    // use local arrays instead of heteregeneous arrays as we want ops on them

    IOStream inFile( fileName, std::ios::in );

    inFile.readFormatted( array, size );

    if ( nEntries != size )
    {
        hmemo::HArray<ValueType> block( nEntries );
        hmemo::ContextPtr ctx = hmemo::Context::getHostPtr();
        SCAI_LOG_DEBUG( logger, "read block first = " << first << ", n = " << nEntries << " from array " << array )

        IndexType inc = 1;
        HArrayUtils::setArraySection( block, 0, inc, array, first, inc, nEntries, common::BinaryOp::COPY, ctx );

        array.swap( block );
    }
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data, outputFileName );
}

void TextIO::readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    readArray( data, inputFileName );
    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "Text does not support multidimensional array, take default shape " << grid )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const std::string& fileName )
{
    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray, fileName, 0, invalidIndex );
    size = denseArray.size();
    ValueType zeroValue = 0;
    HArrayUtils::buildSparseArray( values, indexes, denseArray, zeroValue );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    auto coo = convert<COOStorage<ValueType>>( storage );

    IndexType numRows;
    IndexType numCols;

    HArray<IndexType> cooIA;
    HArray<IndexType> cooJA;
    HArray<ValueType> cooValues;

    coo.splitUp( numRows, numCols, cooIA, cooJA, cooValues );

    IOStream outFile( fileName, std::ios::out );

    int precIndex = 0;
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    if ( mScalarTypeData == common::ScalarType::PATTERN )
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
void TextIO::readData(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>* values,
    const IndexType nnz,
    const std::string& fileName )
{
    HArray<DefaultReal> dIA;
    HArray<DefaultReal> dJA;

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

    HArrayUtils::assign( ia, dIA );  // conversion from DefaultReal to IndexType
    HArrayUtils::assign( ja, dJA );  // conversion from DefaultReal to IndexType

    IndexType minRowIndex = HArrayUtils::reduce( ia, common::BinaryOp::MIN );

    if ( minRowIndex == 0 )
    {
        // okay, seems that indexing start with 0
    }
    else if ( minRowIndex == 1 )
    {
        // offset base = 1, convert it to 0

        HArrayUtils::setScalar( ia, IndexType( 1 ), common::BinaryOp::SUB );
        HArrayUtils::setScalar( ja, IndexType( 1 ), common::BinaryOp::SUB );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ERROR reading file " << fileName << ": minimal row index " << minRowIndex << " is illegal" )
    }
}

/* --------------------------------------------------------------------------------- */

void TextIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName )
{
    IndexType nEntries;

    checkTextFile( numValues, nEntries, fileName.c_str() );

    // As there is no header, we have to read the full file, at least the index values

    HArray<IndexType> ia;
    HArray<IndexType> ja;

    readData<DefaultReal>( ia, ja, NULL, numValues, fileName );

    numRows    = HArrayUtils::max( ia ) + 1;
    numColumns = HArrayUtils::max( ja ) + 1;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName,
    const IndexType firstRow,
    const IndexType nRows )
{
    // binary mode does not matter as we have always formatted output

    // read coo entries lines by line, similiar to _MatrixMarket
    // i , j, val

    IndexType nnz;
    IndexType k;

    checkTextFile( nnz, k, fileName.c_str() );

    SCAI_LOG_INFO( logger, "File : " << fileName << ", #lines = " << nnz << ", #entries = " << k )

    if ( nnz == 0 )
    {
        storage.clear();
        return;
    }

    bool readPattern = mScalarTypeData == common::ScalarType::PATTERN;

    IndexType nEntries = 2;

    if ( !readPattern )
    {
        nEntries = 3;
    }

    SCAI_ASSERT_GE( k, nEntries, "#entries/row in file " << fileName << " must be at least " << nEntries )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    HArray<IndexType> ia;
    HArray<IndexType> ja;
    HArray<ValueType> val;

    if ( readPattern )
    {
        readData<ValueType>( ia, ja, NULL, nnz, fileName );
        val.setSameValue( nnz, ValueType( 1 ) );
    }
    else
    {
        readData<ValueType>( ia, ja, &val, nnz, fileName );
    }

    // we shape the matrix by maximal appearing indexes

    int nrows = HArrayUtils::max( ia ) + 1;
    int ncols = HArrayUtils::max( ja ) + 1;

    COOStorage<ValueType> coo( nrows, ncols, ia, ja, val );

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

void TextIO::writeStorage( const _MatrixStorage& storage, const std::string& fileName )
{
    IOWrapper<TextIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( ( TextIO& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readStorage(
    _MatrixStorage& storage,
    const std::string& fileName,
    const IndexType offsetRow,
    const IndexType nRows )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( ( TextIO& ) *this, storage, fileName, offsetRow, nRows );
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeArray( const hmemo::_HArray& array, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( ( TextIO& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( ( TextIO& ) *this, n, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readArray( hmemo::_HArray& array, const std::string& fileName, const IndexType offset, const IndexType n )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( ( TextIO& ) *this, array, fileName, offset, n );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( ( TextIO& ) *this, size, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

std::string TextIO::getMatrixFileSuffix() const
{
    return TextIO::createValue();
}

/* --------------------------------------------------------------------------------- */

std::string TextIO::getVectorFileSuffix() const
{
    return TextIO::createValue();
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
