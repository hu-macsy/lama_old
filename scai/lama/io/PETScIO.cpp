/*
 * @file PETScIO.cpp
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
 * @brief Implementation of IO methods for PETSc format
 * @author Thomas Brandes
 * @date 10.06.2016
 */


#include <scai/lama/io/PETScIO.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/tracing.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Grid.hpp>

#define PETSC_SUFFIX ".psc"

/** Internal id as specified by PETSc */

#define MAT_FILE_CLASSID 1211216

/** Internal id as specified by PETSc */

#define VEC_FILE_CLASSID 1211214

namespace scai
{

using namespace hmemo;

namespace lama
{

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( PETScIO::logger, "FileIO.PETScIO" )

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* PETScIO::create()
{
    return new PETScIO();
}

std::string PETScIO::createValue()
{
    return PETSC_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool PETScIO::isSupportedMode( const FileMode mode ) const
{
    // only binary is supported


    if ( mode == FORMATTED )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeAt( std::ostream& stream ) const
{
    stream << "PETScIO ( ";
    stream << "suffix = " << PETSC_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only binary, BIG Endian )";
}

/* --------------------------------------------------------------------------------- */

PETScIO::PETScIO()
{
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not available for PETScIO" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    IndexType nrows = array.size();

    std::ios::openmode flags = std::ios::out | std::ios::binary;

    if ( mAppendMode )
    {
        flags |= std::ios::app;
    }
    else
    {
        flags |= std::ios::trunc;
    }

    IOStream outFile( fileName, flags, IOStream::BIG );

    SCAI_LOG_INFO( logger, "File " << fileName << " now open for binary write, append = " << mAppendMode )

    HArray<IndexType> headValues( { VEC_FILE_CLASSID, nrows } );

    outFile.writeBinary( headValues, mScalarTypeIndex );
    outFile.writeBinary( array, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeSparseImpl(
    const IndexType size,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values,
    const std::string& fileName )
{
    // sparse unsupported for this file format, write it dense

    HArray<ValueType> denseArray;
    ValueType zero = 0;
    utilskernel::HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );
    writeArrayImpl( denseArray, fileName );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readArrayInfo( IndexType& size, const std::string& fileName )
{
    SCAI_REGION( "IO.PETSc.readArrayInfo" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    SCAI_LOG_INFO( logger, "Read array info from file " << fileName )

    IOStream inFile( fileName, flags, IOStream::BIG );

    HArray<IndexType> headerVals;

    inFile.readBinary( headerVals, 2, common::ScalarType::INDEX_TYPE );

    IndexType classid = headerVals[0];

    SCAI_ASSERT_EQUAL( VEC_FILE_CLASSID, classid, "illegal VEC_FILE_CLASSID" )

    size = headerVals[1];
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName,
    const IndexType first,
    const IndexType n )
{
    SCAI_REGION( "IO.PETSc.readArray" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    SCAI_LOG_INFO( logger, "Read array<" << common::TypeTraits<ValueType>::id()
                   << "> from file " << fileName << ", type = " << mScalarTypeData )

    IOStream inFile( fileName, flags, IOStream::BIG );

    HArray<IndexType> headerVals;

    inFile.readBinary( headerVals, 2, common::ScalarType::INDEX_TYPE );

    IndexType classid = headerVals[0];
    IndexType size    = headerVals[1];

    SCAI_LOG_INFO( logger, "Read: id = " << classid << ", size = " << size )

    SCAI_ASSERT_EQUAL( VEC_FILE_CLASSID, classid, "illegal VEC_FILE_CLASSID" )

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

    // check if the specified data size fits the expected data type

    inFile.readBinary( array, size, mScalarTypeData );

    inFile.close();

    if ( nEntries != array.size() )
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
void PETScIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const std::string& fileName )
{
    SCAI_REGION( "IO.PETSc.readSparse" )

    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray, fileName, 0, invalidIndex );
    size = denseArray.size();
    ValueType zero = 0;
    utilskernel::HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_REGION( "IO.PETSc.writeStorage" )

    SCAI_ASSERT( mFileMode != FORMATTED, "Formatted output not available for PETScIO" )

    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    IndexType nrows = storage.getNumRows();
    IndexType ncols = storage.getNumColumns();

    HArray<IndexType> csrIA;    // first offsets, later sizes
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    IndexType nnz = csrJA.size();

    // we need the CSR sizes, not the offsets

    utilskernel::HArrayUtils::unscan( csrIA );

    std::ios::openmode flags = std::ios::out | std::ios::binary;

    if ( mAppendMode )
    {
        flags |= std::ios::app;
    }
    else
    {
        flags |= std::ios::trunc;
    }

    IOStream outFile( fileName, flags, IOStream::BIG );

    SCAI_LOG_INFO( logger, "File " << fileName << " now open for binary write, append = " << mAppendMode )

    // Note: PETSc starts indexing with 0

    HArray<IndexType> headValues( 4 );

    headValues[0] = MAT_FILE_CLASSID;
    headValues[1] = nrows;
    headValues[2] = ncols;
    headValues[3] = nnz;

    // for binary output we make conversions to mScalarTypeData, mScalarTypeIndex

    outFile.writeBinary( headValues, mScalarTypeIndex );
    outFile.writeBinary( csrIA, mScalarTypeIndex );
    outFile.writeBinary( csrJA , mScalarTypeIndex );

    // output of values is skipped for PATTERN

    if ( mScalarTypeData != common::ScalarType::PATTERN )
    {
        outFile.writeBinary( csrValues, mScalarTypeData );
    }
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName )
{
    SCAI_REGION( "IO.PETSc.readStorageInfo" )

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    SCAI_LOG_INFO( logger, "Read storage info from file " << fileName )

    IOStream inFile( fileName, flags, IOStream::BIG );

    HArray<IndexType> headerVals;

    inFile.readBinary( headerVals, 4, common::TypeTraits<IndexType>::stype );

    IndexType classid = headerVals[0];

    SCAI_ASSERT_EQUAL( MAT_FILE_CLASSID, classid, "illegal MAT_FILE_CLASSID" )

    numRows      = headerVals[1];
    numColumns   = headerVals[2];
    numValues    = headerVals[3];
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName,
    const IndexType firstRow,
    const IndexType nRows )
{
    SCAI_REGION( "IO.PETSc.readStorage" )

    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    SCAI_LOG_INFO( logger, "Read storage<" << common::TypeTraits<ValueType>::id() << "> from file " << fileName )

    IOStream inFile( fileName, flags, IOStream::BIG );

    HArray<IndexType> headerVals;

    inFile.readBinary( headerVals, 4, common::TypeTraits<IndexType>::stype );

    IndexType classid = headerVals[0];
    IndexType numRows   = headerVals[1];
    IndexType numCols   = headerVals[2];
    IndexType nnz     = headerVals[3];

    SCAI_LOG_INFO( logger, "Read: id = " << MAT_FILE_CLASSID << ", #rows = " << numRows
                   << ", #cols = " << numCols << ", #nnz = " << nnz )

    SCAI_ASSERT_EQUAL( MAT_FILE_CLASSID, classid, "illegal MAT_FILE_CLASSID" )

    HArray<IndexType> csrSizes;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    inFile.readBinary( csrSizes, numRows, mScalarTypeIndex );
    inFile.readBinary( csrJA, nnz, mScalarTypeIndex );

    if ( mScalarTypeData != common::ScalarType::PATTERN )
    {
        inFile.readBinary( csrValues, nnz, mScalarTypeData );
    }
    else
    {
        csrValues.setSameValue( nnz, ValueType( 1 ) );
    }

    if ( firstRow == 0 && nRows == invalidIndex )
    {
        storage.setCSRData( numRows, numCols, csrSizes, csrJA, csrValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "read block not yet available" )
    }
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data, outputFileName );
}

void PETScIO::readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& inputFileName )
{
    readArray( data, inputFileName );
    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "PETSc does not support multidimensional array, take default shape " << grid )
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeStorage( const _MatrixStorage& storage, const std::string& fileName )
{
    IOWrapper<PETScIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( ( PETScIO& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readStorage(
    _MatrixStorage& storage,
    const std::string& fileName,
    const IndexType offsetRow,
    const IndexType nRows )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( ( PETScIO& ) *this, storage, fileName, offsetRow, nRows );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeArray( const hmemo::_HArray& array, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( ( PETScIO& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( ( PETScIO& ) *this, n, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readArray( hmemo::_HArray& array, const std::string& fileName, const IndexType offset, const IndexType n )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( ( PETScIO& ) *this, array, fileName, offset, n );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const std::string& fileName )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( ( PETScIO& ) *this, size, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

std::string PETScIO::getMatrixFileSuffix() const
{
    return PETScIO::createValue();
}

/* --------------------------------------------------------------------------------- */

std::string PETScIO::getVectorFileSuffix() const
{
    return PETScIO::createValue();
}
/* --------------------------------------------------------------------------------- */

#define SCAI_PETSC_METHOD_INSTANTIATIONS( _type )           \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::writeArrayImpl(                           \
        const hmemo::HArray<_type>& array,                  \
        const std::string& fileName );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::readArrayImpl(                            \
        hmemo::HArray<_type>& array,                        \
        const std::string& arrayFileName,                   \
        const IndexType ,                                   \
        const IndexType );                                  \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::writeSparseImpl(                          \
        const IndexType size,                               \
        const HArray<IndexType>& index,                     \
        const HArray<_type>& values,                        \
        const std::string& fileName );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::readSparseImpl(                           \
        IndexType& size,                                    \
        HArray<IndexType>& indexes,                         \
        HArray<_type>& values,                              \
        const std::string& fileName );         

SCAI_COMMON_LOOP( SCAI_PETSC_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_PETSC_METHOD_INSTANTIATIONS

#define SCAI_PETSC_METHOD_INSTANTIATIONS( _type )       \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void PETScIO::writeStorageImpl(                     \
        const MatrixStorage<_type>& storage,            \
        const std::string& fileName );                  \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void PETScIO::readStorageImpl(                      \
        MatrixStorage<_type>& storage,                  \
        const std::string& matrixFileName,              \
        const IndexType firstRow,                       \
        const IndexType nRows );                     

SCAI_COMMON_LOOP( SCAI_PETSC_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_PETSC_METHOD_INSTANTIATIONS

}  // lama

}  // scai
