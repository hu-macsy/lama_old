/**
 * @file TextIO.cpp
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

    if ( mode == FileMode::BINARY )
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
void TextIO::checkTextFile( IndexType& nLines, IndexType& nEntries )
{
    SCAI_LOG_DEBUG( logger, "checkTextFile " << mFile.getFileName() )

    nLines   = 0;
    nEntries = 0;

    std::string line;
    std::vector<std::string> tokens;

    while ( std::getline( mFile, line ) )
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

    SCAI_LOG_INFO( logger, "checkTextFile " << mFile.getFileName() << ": #lines = " << nLines << ", #entries = " << nEntries )

    mFile.clear();
    mFile.seekg( 0, std::ios::beg );
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( TextIO::logger, "FileIO.TextIO" )

/* --------------------------------------------------------------------------------- */

void TextIO::openIt( const std::string& fileName, const char* openMode )
{
    if ( strcmp( openMode, "w" ) == 0 )
    {
        SCAI_ASSERT( mFileMode != FileMode::BINARY, "Binary mode not supported for " << *this )
        mFile.open( fileName, std::ios::out | std::ios::trunc );
    }
    else if ( strcmp( openMode, "r" ) == 0 )
    {
        mFile.open( fileName, std::ios::in );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported open mode for Text file: " << openMode )
    }
}

/* --------------------------------------------------------------------------------- */

void TextIO::closeIt()
{
    mFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeArrayImpl( const hmemo::HArray<ValueType>& array )
{
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    mFile.writeFormatted( array, precData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeSparseImpl(
    const IndexType size,
    const ValueType& zero,
    const hmemo::HArray<IndexType>& indexes,
    const hmemo::HArray<ValueType>& values )
{
    // Note: size is ignored for TextIO of sparse vector

    HArray<ValueType> denseArray;
    HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );
    writeArrayImpl( denseArray );

    // Better solution to write sparse data, but not unique when reading, no header
    //     mFile.writeFormatted( indexes, precIndexes, values, precData );
}

/* --------------------------------------------------------------------------------- */

void TextIO::getArrayInfo( IndexType& size )
{
    IndexType nEntries;   // dummy variable needed for checkTextFile

    // each array entry in one line, so count the number of lines in the file

    checkTextFile( size, nEntries );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readArrayImpl( hmemo::HArray<ValueType>& array )
{
    IndexType size;   // number of lines, size of array
    IndexType k;      // number of entries in one line

    checkTextFile( size, k );

    SCAI_ASSERT_LE( k, 2, "#entries/row in file " << mFile.getFileName() << " must not excced 2" )

    mFile.readFormatted( array, size );
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data );
}

void TextIO::readGridArray( hmemo::_HArray& data, common::Grid& grid )
{
    readArray( data );
    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "Text does not support multidimensional array, take default shape " << grid )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readSparseImpl(
    IndexType& size,
    ValueType& zero,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values )
{
    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray );
    size = denseArray.size();
    zero = 0;
    HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::writeStorageImpl( const MatrixStorage<ValueType>& storage )
{
    auto coo = convert<COOStorage<ValueType>>( storage );

    IndexType numRows;
    IndexType numCols;

    HArray<IndexType> cooIA;
    HArray<IndexType> cooJA;
    HArray<ValueType> cooValues;

    coo.splitUp( numRows, numCols, cooIA, cooJA, cooValues );

    int precIndex = 0;
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    if ( mScalarTypeData == common::ScalarType::PATTERN )
    {
        mFile.writeFormatted( cooIA, precIndex, cooJA, precIndex );
    }
    else
    {
        mFile.writeFormatted( cooIA, precIndex, cooJA, precIndex, cooValues, precData );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readData(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>* values,
    const IndexType nnz )
{
    HArray<DefaultReal> dIA;
    HArray<DefaultReal> dJA;

    if ( values != NULL )
    {
        mFile.readFormatted( dIA, dJA, *values, nnz );
    }
    else
    {
        mFile.readFormatted( dIA, dJA, nnz );
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
        COMMON_THROWEXCEPTION( "ERROR reading file " << mFile.getFileName() << ": minimal row index " << minRowIndex << " is illegal" )
    }
}

/* --------------------------------------------------------------------------------- */

void TextIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    IndexType nEntries;

    checkTextFile( numValues, nEntries );

    // As there is no header, we have to read the full file, at least the index values

    HArray<IndexType> ia;
    HArray<IndexType> ja;

    readData<DefaultReal>( ia, ja, NULL, numValues );

    numRows    = HArrayUtils::max( ia ) + 1;
    numColumns = HArrayUtils::max( ja ) + 1;

    mFile.clear();
    mFile.seekg( 0, std::ios::beg );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void TextIO::readStorageImpl( MatrixStorage<ValueType>& storage )
{
    // binary mode does not matter as we have always formatted output

    // read coo entries lines by line, similiar to _MatrixMarket
    // i , j, val

    IndexType nnz;
    IndexType k;

    checkTextFile( nnz, k );

    SCAI_LOG_INFO( logger, "File : " << mFile.getFileName() << ", #lines = " << nnz << ", #entries = " << k )

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

    SCAI_ASSERT_GE( k, nEntries, "#entries/row in file " << mFile.getFileName() << " must be at least " << nEntries )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    HArray<IndexType> ia;
    HArray<IndexType> ja;
    HArray<ValueType> val;

    if ( readPattern )
    {
        readData<ValueType>( ia, ja, NULL, nnz );
        val.setSameValue( nnz, ValueType( 1 ) );
    }
    else
    {
        readData<ValueType>( ia, ja, &val, nnz );
    }

    // we shape the matrix by maximal appearing indexes

    int nrows = HArrayUtils::max( ia ) + 1;
    int ncols = HArrayUtils::max( ja ) + 1;

    COOStorage<ValueType> coo( nrows, ncols, ia, ja, val );

    storage = coo;
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeStorage( const _MatrixStorage& storage )
{
    IOWrapper<TextIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readStorage( _MatrixStorage& storage )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeArray( const hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

void TextIO::writeSparse( const IndexType n, const void* zero, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparse( *this, n, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readArray( hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

void TextIO::readSparse( IndexType& size, void* zero, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<TextIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparse( *this, size, zero, indexes, values );
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
