/**
 * @file PETScIO.cpp
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

    if ( mode == FileMode::FORMATTED )
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

void PETScIO::open( const char* fileName, const char* fileMode )
{
    SCAI_ASSERT( mFileMode != FileMode::FORMATTED, "Formatted output not available for PETScIO" )

    std::ios::openmode flags = std::ios::binary;

    if ( strcmp( fileMode, "w" ) == 0 )
    {
        flags |= std::ios::out | std::ios::trunc;
    }
    else if ( strcmp( fileMode, "a" ) == 0 )
    {
        flags |= std::ios::out | std::ios::app;
    }
    else if ( strcmp( fileMode, "r" ) == 0 )
    {
        flags |= std::ios::in ;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported file mode for PETSc file: " << fileMode )
    }

    mFile.open( fileName, flags, IOStream::BIG );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::close()
{
    mFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeArrayImpl( const hmemo::HArray<ValueType>& array )
{
    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    IndexType nrows = array.size();

    HArray<IndexType> headValues( { VEC_FILE_CLASSID, nrows } );

    mFile.writeBinary( headValues, mScalarTypeIndex );
    mFile.writeBinary( array, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeSparseImpl(
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

void PETScIO::getArrayInfo( IndexType& size )
{
    SCAI_REGION( "IO.PETSc.getArrayInfo" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    SCAI_LOG_INFO( logger, "Read array info from file " << mFile.getFileName() )

    HArray<IndexType> headerVals;

    std::streampos pos = mFile.tellg();

    mFile.readBinary( headerVals, 2, common::ScalarType::INDEX_TYPE );

    mFile.clear();       // important to reset flags
    mFile.seekg( pos );

    IndexType classid = headerVals[0];

    SCAI_ASSERT_EQUAL( VEC_FILE_CLASSID, classid, "illegal VEC_FILE_CLASSID" )

    size = headerVals[1];
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::readArrayImpl( hmemo::HArray<ValueType>& array )
{
    SCAI_REGION( "IO.PETSc.readArray" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    SCAI_LOG_INFO( logger, "Read array<" << common::TypeTraits<ValueType>::id()
                   << "> from file " << mFile.getFileName() << ", type = " << mScalarTypeData )

    HArray<IndexType> headerVals;

    mFile.readBinary( headerVals, 2, common::ScalarType::INDEX_TYPE );

    IndexType classid = headerVals[0];
    IndexType size    = headerVals[1];

    SCAI_LOG_INFO( logger, "Read: id = " << classid << ", size = " << size )

    SCAI_ASSERT_EQUAL( VEC_FILE_CLASSID, classid, "illegal VEC_FILE_CLASSID" )

    // check if the specified data size fits the expected data type

    mFile.readBinary( array, size, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::readSparseImpl(
    IndexType& size,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values )
{
    SCAI_REGION( "IO.PETSc.readSparse" )

    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray );
    size = denseArray.size();
    ValueType zero = 0;
    utilskernel::HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::writeStorageImpl( const MatrixStorage<ValueType>& storage )
{
    SCAI_REGION( "IO.PETSc.writeStorage" )

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

    // Note: PETSc starts indexing with 0

    HArray<IndexType> headValues( 4 );

    headValues[0] = MAT_FILE_CLASSID;
    headValues[1] = nrows;
    headValues[2] = ncols;
    headValues[3] = nnz;

    // for binary output we make conversions to mScalarTypeData, mScalarTypeIndex

    mFile.writeBinary( headValues, mScalarTypeIndex );
    mFile.writeBinary( csrIA, mScalarTypeIndex );
    mFile.writeBinary( csrJA , mScalarTypeIndex );

    // output of values is skipped for PATTERN

    if ( mScalarTypeData != common::ScalarType::PATTERN )
    {
        mFile.writeBinary( csrValues, mScalarTypeData );
    }
}

/* --------------------------------------------------------------------------------- */

void PETScIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    SCAI_REGION( "IO.PETSc.getStorageInfo" )

    HArray<IndexType> headerVals;

    std::streampos pos = mFile.tellg();

    SCAI_LOG_INFO( logger, "get storage info from file " << mFile.getFileName() << ", pos = " << pos )

    mFile.readBinary( headerVals, 4, common::TypeTraits<IndexType>::stype );

    mFile.clear();       // important to reset flags
    mFile.seekg( pos );

    IndexType classid = headerVals[0];

    SCAI_ASSERT_EQUAL( MAT_FILE_CLASSID, classid, "illegal MAT_FILE_CLASSID" )

    numRows      = headerVals[1];
    numColumns   = headerVals[2];
    numValues    = headerVals[3];
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PETScIO::readStorageImpl( MatrixStorage<ValueType>& storage )
{
    SCAI_REGION( "IO.PETSc.readStorage" )

    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    SCAI_LOG_INFO( logger, "Read storage<" << common::TypeTraits<ValueType>::id() << "> from file " << mFile.getFileName() )

    HArray<IndexType> headerVals;

    mFile.readBinary( headerVals, 4, common::TypeTraits<IndexType>::stype );

    IndexType classid = headerVals[0];
    IndexType numRows = headerVals[1];
    IndexType numCols = headerVals[2];
    IndexType nnz     = headerVals[3];

    SCAI_LOG_INFO( logger, "Read: id = " << MAT_FILE_CLASSID << ", #rows = " << numRows
                   << ", #cols = " << numCols << ", #nnz = " << nnz )

    SCAI_ASSERT_EQUAL( MAT_FILE_CLASSID, classid, "illegal MAT_FILE_CLASSID" )

    HArray<IndexType> csrSizes;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    mFile.readBinary( csrSizes, numRows, mScalarTypeIndex );
    mFile.readBinary( csrJA, nnz, mScalarTypeIndex );

    if ( mScalarTypeData != common::ScalarType::PATTERN )
    {
        mFile.readBinary( csrValues, nnz, mScalarTypeData );
    }
    else
    {
        csrValues.setSameValue( nnz, ValueType( 1 ) );
    }

    storage.setCSRData( numRows, numCols, csrSizes, csrJA, csrValues );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data );
}

void PETScIO::readGridArray( hmemo::_HArray& data, common::Grid& grid )
{
    readArray( data );
    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "PETSc does not support multidimensional array, take default shape " << grid )
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeStorage( const _MatrixStorage& storage )
{
    IOWrapper<PETScIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readStorage( _MatrixStorage& storage )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeArray( const hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( *this, array );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( *this, n, indexes, values );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readArray( hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( *this, array );
}

/* --------------------------------------------------------------------------------- */

void PETScIO::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<PETScIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( *this, size, indexes, values );
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
        const hmemo::HArray<_type>& array );                \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::readArrayImpl(                            \
        hmemo::HArray<_type>& array );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::writeSparseImpl(                          \
        const IndexType size,                               \
        const HArray<IndexType>& index,                     \
        const HArray<_type>& values );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void PETScIO::readSparseImpl(                           \
        IndexType& size,                                    \
        HArray<IndexType>& indexes,                         \
        HArray<_type>& values );

SCAI_COMMON_LOOP( SCAI_PETSC_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_PETSC_METHOD_INSTANTIATIONS

#define SCAI_PETSC_METHOD_INSTANTIATIONS( _type )       \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void PETScIO::writeStorageImpl(                     \
        const MatrixStorage<_type>& storage );          \
                                                        \
    template COMMON_DLL_IMPORTEXPORT                    \
    void PETScIO::readStorageImpl(                      \
        MatrixStorage<_type>& storage );                \

SCAI_COMMON_LOOP( SCAI_PETSC_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_PETSC_METHOD_INSTANTIATIONS

}  // lama

}  // scai
