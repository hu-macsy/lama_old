/**
 * @file LamaIO.cpp
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
 * @brief Implementation of IO methods for Lama format
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#include <scai/lama/io/LamaIO.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/io/IOWrapper.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/tracing.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Grid.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <iomanip>

#define LAMA_SUFFIX ".lmf"

namespace scai
{

using common::ScalarType;
using common::Grid;

using namespace hmemo;

using utilskernel::HArrayUtils;

namespace lama
{

std::ostream& operator<<( std::ostream& stream, const LamaClassId& object )
{
    switch ( object )
    {
        case LamaClassId::SPARSE_VECTOR:
            stream << "SPARSE_VECTOR";
            break;

        case LamaClassId::GRID_VECTOR:
            stream << "GRID_VECTOR";
            break;

        case LamaClassId::CSR_MATRIX:
            stream << "CSR_MATRIX";
            break;

        case LamaClassId::COO_MATRIX:
            stream << "COO_MATRIX";
            break;

        default:
            stream << static_cast<int>( object );
    }

    stream << " ( " << std::hex << static_cast<int>( object ) << std::dec << " )";

    return stream;
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( LamaIO::logger, "FileIO.LamaIO" )

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* LamaIO::create()
{
    return new LamaIO();
}

std::string LamaIO::createValue()
{
    return LAMA_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool LamaIO::isSupportedMode( const FileMode mode ) const
{
    // only binary is supported

    if ( mode == FileMode::FORMATTED )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void LamaIO::writeAt( std::ostream& stream ) const
{
    stream << "LamaIO ( ";
    stream << "suffix = " << LAMA_SUFFIX << ", ";
    stream << "name = " << mFileName << ", ";
    FileIO::writeMode( stream );
    stream << ", only binary )";
}

/* --------------------------------------------------------------------------------- */

LamaIO::LamaIO()
{
    SCAI_LOG_INFO( logger, "Constructor: " << *this )
}

/* --------------------------------------------------------------------------------- */

void LamaIO::openIt( const std::string& fileName, const char* fileMode )
{
    SCAI_ASSERT( mFileMode != FileMode::FORMATTED, "Formatted output not available for LamaIO" )

    auto comm = getDistributedIOMode()  == DistributedIOMode::COLLECTIVE
                ? dmemo::Communicator::getCommunicatorPtr()
                : dmemo::Communicator::getCommunicatorPtr( dmemo::Communicator::NO );

    mFile = comm->collectiveFile();

    SCAI_ASSERT_ERROR( mFile.get(), "Could not get collective file object for comm = " << comm )

    mFileName = fileName;

    mFile->open( fileName.c_str(), fileMode );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::closeIt()
{
    mFile->close();

    mFile.reset();
    mFileName.clear();
}

/* --------------------------------------------------------------------------------- */

bool LamaIO::hasCollectiveIO() const
{
    return true;
}

/* --------------------------------------------------------------------------------- */

void LamaIO::writeHeader( const LamaClassId classId )
{
    const int header[] = { static_cast<int>( classId ),
                           static_cast<int>( mScalarTypeIndex ),
                           static_cast<int>( mScalarTypeData )
                         };

    mFile->writeSingle( header, 3 );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeArrayImpl( const hmemo::HArray<ValueType>& array )
{
    // handle array like a one-dimensional grid array

    common::Grid1D grid1d( array.size() );
    writeGridImpl( array, grid1d );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeSparseImpl(
    const IndexType size,
    const ValueType& zero,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values )
{
    SCAI_ASSERT_EQ_ERROR( indexes.size(), values.size(), "mismatch cooridinates/values sizes" )

    // int       SPARSE_VECTOR
    // int       IndexType.id
    // int       DataType.id
    // IndexType vector size
    // DataType  zero value
    // IndexType number of non-zeros
    // IndexType non-zero indexes
    // DataType  non-zero values

    const auto& comm = mFile->getCommunicator();

    IndexType n   = comm.sum( size );
    IndexType nnz = comm.sum( indexes.size() );

    writeHeader( LamaClassId::SPARSE_VECTOR );

    mFile->writeSingle( n, mScalarTypeIndex ),
          mFile->writeSingle( zero, mScalarTypeData );
    mFile->writeSingle( nnz, mScalarTypeIndex ),
          mFile->writeAll( indexes, mScalarTypeIndex ),
          mFile->writeAll( values, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::getArrayInfo( IndexType& size )
{
    SCAI_REGION( "IO.Lama.getArrayInfo" )

    SCAI_LOG_INFO( logger, "Read array info from file " << mFileName )

    size_t pos = mFile->getOffset();

    int header[3];     // array to read the first three entries from the file

    mFile->readSingle( header, 3 );

    LamaClassId classId = LamaClassId( header[0] );

    auto fileIndexType = ScalarType( header[1] );

    if ( classId == LamaClassId::SPARSE_VECTOR )
    {
        mFile->readSingle( size, fileIndexType );
    }
    else if ( classId == LamaClassId::GRID_VECTOR )
    {
        size = 1;

        IndexType nDims;
        mFile->readSingle( nDims, fileIndexType );
        for ( IndexType i = 0; i < nDims; ++i )
        {
            IndexType dimSize;
            mFile->readSingle( dimSize, fileIndexType );
            size *= dimSize;
        }
    }
    else
    {
        size = invalidIndex;
    }

    mFile->setOffset( pos );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readArrayImpl( hmemo::HArray<ValueType>& array )
{
    SCAI_REGION( "IO.Lama.readArray" )

    common::Grid grid;

    readGridImpl( array, grid );

    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "multidimensional array " << grid << ") in file " 
                                << *this << " considered as 1D" )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readSparseImpl(
    IndexType& size,
    ValueType& zero,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values )
{
    // int       SPARSE_VECTOR
    // int       IndexType.id
    // int       DataType.id
    // IndexType vector size
    // DataType  zero value
    // IndexType number of non-zeros
    // IndexType non-zero indexes
    // DataType  non-zero values

    const auto& comm = mFile->getCommunicator();

    int header[3];     // array to read the header of next entry

    mFile->readSingle( header, 3 );

    SCAI_ASSERT_EQ_ERROR( LamaClassId( header[0] ), LamaClassId::SPARSE_VECTOR, "no sparse vector in file" )

    auto fileIndexType = ScalarType( header[1] );
    auto fileDataType  = ScalarType( header[2] );

    IndexType numValues;   // total number of non-zeros

    mFile->readSingle( size, fileIndexType );
    mFile->readSingle( zero, fileDataType );
    mFile->readSingle( numValues, fileIndexType );

    // each processor assembles a chunk of entries

    auto dist = dmemo::blockDistribution( numValues, mFile->getCommunicatorPtr() );
    auto localNumValues = dist->getLocalSize();

    mFile->readAll( indexes, localNumValues, dist->lb(), fileIndexType );
    mFile->readAll( values, localNumValues, dist->lb(), fileDataType );

    IndexType allNNZ = comm.sum( localNumValues );
    SCAI_ASSERT_EQ_ERROR( allNNZ, numValues, "serious mismatch for parallel read CSR matrix data" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeDense( const DenseStorage<ValueType>& dense )
{
    // dense storage is written as 2D grid data

    common::Grid2D grid( dense.getNumRows(), dense.getNumColumns() );
    writeGridImpl( dense.getValues(), grid );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readDense( DenseStorage<ValueType>& dense )
{
    // dense storage is written as 2D grid data

    HArray<ValueType> data;
    common::Grid grid;
    IndexType dummyRows;
    IndexType dummyCols;

    dense.splitUp( dummyRows, dummyCols, data );

    readGridImpl( data, grid );

    SCAI_ASSERT_EQ_ERROR( grid.nDims(), 2, "grid data is not general matrix" )

    dense = DenseStorage<ValueType>( grid.size( 0 ), grid.size( 1 ), std::move( data ) );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeCSR( const CSRStorage<ValueType>& csr )
{
    const auto& comm = mFile->getCommunicator();

    const IndexType numRows = comm.sum( csr.getNumRows() );
    const IndexType numCols = csr.getNumColumns();
    const IndexType numValues = comm.sum( csr.getNumValues() );

    // Meta-data for CSR storage, header + 'global' sizes

    writeHeader( LamaClassId::CSR_MATRIX );

    mFile->writeSingle( numRows, mScalarTypeIndex );
    mFile->writeSingle( numCols, mScalarTypeIndex );
    mFile->writeSingle( numValues, mScalarTypeIndex );

    HArray<IndexType> csrSizes = csr.getIA();

    HArrayUtils::unscan( csrSizes );  // we write sizes into file, not offsets

    mFile->writeAll( csrSizes, mScalarTypeIndex );
    mFile->writeAll( csr.getJA(), mScalarTypeIndex );

    if ( mScalarTypeData != ScalarType::PATTERN )
    {
        mFile->writeAll( csr.getValues(), mScalarTypeData );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readCSR( CSRStorage<ValueType>& csr )
{
    const auto& comm = mFile->getCommunicator();

    int val[3];     // array to read the first three entries from the file

    mFile->readSingle( val, 3 );

    SCAI_LOG_INFO( logger, comm << ": read header for CSR matrix: " << val[0] << ", " << val[1] << ", " << val[2] )

    SCAI_ASSERT_EQ_ERROR( LamaClassId( val[0] ), LamaClassId::CSR_MATRIX, "no CSR matrix in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType numRows;
    IndexType numCols;
    IndexType numValues;

    mFile->readSingle( numRows, fileIndexType );
    mFile->readSingle( numCols, fileIndexType );
    mFile->readSingle( numValues, fileIndexType );

    SCAI_LOG_DEBUG( logger, comm << ": read CSR matrix " << numRows << " x " << numCols << ", NNZ = " << numValues )

    auto dist = dmemo::blockDistribution( numRows, mFile->getCommunicatorPtr() );
    auto localNumRows = dist->getLocalSize();

    HArray<IndexType> csrIA( localNumRows + 1 );

    mFile->readAll( csrIA, localNumRows, dist->lb(), fileIndexType );

    const IndexType localNNZ = HArrayUtils::scan1( csrIA );

    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    mFile->readAll( csrJA, localNNZ, fileIndexType );

    if ( fileDataType != ScalarType::PATTERN )
    {
        mFile->readAll( csrValues, localNNZ, fileDataType );
    }
    else if ( mScalarTypeData == common::ScalarType::PATTERN )
    {
        csrValues.setSameValue( localNNZ, ValueType( 1 ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "pattern CSR in file, no data, set SCAI_IO_DATA_TYPE=PATTERN" )
    }

    csr = CSRStorage<ValueType>( localNumRows, numCols, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );

    IndexType allNNZ = comm.sum( localNNZ );
    SCAI_ASSERT_EQ_ERROR( allNNZ, numValues, "serious mismatch for parallel read CSR matrix data" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeStorageImpl( const MatrixStorage<ValueType>& storage )
{
    SCAI_REGION( "IO.Lama.writeStorage" )

    if ( storage.getFormat() == Format::DENSE )
    {
        writeDense( static_cast<const DenseStorage<ValueType>&>( storage ) );
    }
    else if ( storage.getFormat() == Format::CSR )
    {
        writeCSR( static_cast<const CSRStorage<ValueType>&>( storage ) );
    }
    else
    {
        CSRStorage<ValueType> tmp;
        tmp.assign( storage );
        writeCSR( tmp );
    }
}

/* --------------------------------------------------------------------------------- */

LamaClassId LamaIO::getClassId()
{
    size_t pos = mFile->getOffset();

    int classId;

    mFile->readSingle( classId );

    mFile->setOffset( pos );

    return LamaClassId( classId );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readStorageImpl( MatrixStorage<ValueType>& storage )
{
    SCAI_REGION( "IO.Lama.readStorage" )

    auto id = getClassId();

    if ( id == LamaClassId::CSR_MATRIX )
    {
        if ( storage.getFormat() == Format::CSR )
        {
            readCSR( static_cast<CSRStorage<ValueType>&>( storage ) );
        }
        else
        {
            CSRStorage<ValueType> tmp;
            readCSR( tmp );
            storage = tmp;
        }
    }
    else if ( id == LamaClassId::GRID_VECTOR )
    {
        if ( storage.getFormat() == Format::DENSE )
        {
            readDense( static_cast<DenseStorage<ValueType>&>( storage ) );
        }
        else
        {
            DenseStorage<ValueType> tmp;
            readStorage( tmp );
            storage = tmp;
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "read storage, but illegal class id: " << id )
    }
}

/* --------------------------------------------------------------------------------- */

void LamaIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    SCAI_REGION( "IO.Lama.getStorageInfo" )

    size_t pos = mFile->getOffset();

    int header[3];   // header data for each entry

    mFile->readSingle( header, 3 );

    SCAI_ASSERT_EQ_ERROR( LamaClassId( header[0] ), LamaClassId::CSR_MATRIX, "no CSR matrix in file" )

    auto fileIndexType = ScalarType( header[1] );

    mFile->readSingle( numRows, fileIndexType );
    mFile->readSingle( numColumns, fileIndexType );
    mFile->readSingle( numValues, fileIndexType );

    mFile->setOffset( pos );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::writeGridImpl( const hmemo::HArray<ValueType>& data, const common::Grid& grid )
{
    writeHeader( LamaClassId::GRID_VECTOR );

    IndexType nDims = grid.nDims();

    const auto& comm = mFile->getCommunicator();

    common::Grid globalGrid( grid );
    globalGrid.setSize( 0, comm.sum( globalGrid.size( 0 ) ) );

    mFile->writeSingle( nDims, mScalarTypeIndex );
    mFile->writeSingle( globalGrid.sizes(), nDims, mScalarTypeIndex );

    mFile->writeAll( data, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void LamaIO::readGridImpl( hmemo::HArray<ValueType>& data, common::Grid& grid )
{
    int val[3];     // array to read the first three entries from the file

    mFile->readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( LamaClassId( val[0] ), LamaClassId::GRID_VECTOR, "no grid vector in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType nDims;

    mFile->readSingle( nDims, fileIndexType );

    SCAI_ASSERT_LE_ERROR( nDims, SCAI_GRID_MAX_DIMENSION, "num dimensions too high" )

    IndexType dims[SCAI_GRID_MAX_DIMENSION];

    mFile->readSingle( dims, nDims, fileIndexType );

    // in case of parallel I/O we assume block distribution in first dimension

    auto dist = dmemo::blockDistribution( dims[0], mFile->getCommunicatorPtr() );

    IndexType ntail = 1;

    for ( IndexType i = 1; i < nDims; ++i )
    {
        ntail *= dims[i];
    }

    dims[0] = dist->getLocalSize();

    grid = Grid( nDims, dims );

    mFile->readAll( data, grid.size(), dist->lb() * ntail, fileDataType );
}

/* -------------------------------------------------------------------------------- */

void LamaIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid )
{
    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeGrid( *this, data, grid );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::readGridArray( hmemo::_HArray& data, common::Grid& grid )
{
    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::readGrid( *this, data, grid );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::writeStorage( const _MatrixStorage& storage )
{
    IOWrapper<LamaIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::readStorage( _MatrixStorage& storage )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<LamaIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::writeArray( const hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::writeSparse( const IndexType n, const void* zero, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparse( *this, n, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::readArray( hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

void LamaIO::readSparse( IndexType& size, void* zero, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<LamaIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparse( *this, size, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

std::string LamaIO::getMatrixFileSuffix() const
{
    return LamaIO::createValue();
}

/* --------------------------------------------------------------------------------- */

std::string LamaIO::getVectorFileSuffix() const
{
    return LamaIO::createValue();
}
/* --------------------------------------------------------------------------------- */

#define SCAI_LAMA_METHOD_INSTANTIATIONS( _type )            \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::writeArrayImpl(                            \
            const HArray<_type>& );                             \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::readArrayImpl(                             \
            HArray<_type>& );                                   \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::writeSparseImpl(                           \
            const IndexType size,                               \
            const _type& zero,                                  \
            const HArray<IndexType>& index,                     \
            const HArray<_type>& values );                      \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::readSparseImpl(                            \
            IndexType& size,                                    \
            _type& zero,                                        \
            HArray<IndexType>& indexes,                         \
            HArray<_type>& values );                            \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::writeGridImpl(                             \
            const HArray<_type>&,                               \
            const common::Grid& );                              \
    \
    template COMMON_DLL_IMPORTEXPORT                        \
    void LamaIO::readGridImpl(                              \
            HArray<_type>& array,                               \
            common::Grid& );

SCAI_COMMON_LOOP( SCAI_LAMA_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_LAMA_METHOD_INSTANTIATIONS

#define SCAI_LAMA_METHOD_INSTANTIATIONS( _type )        \
    \
    template COMMON_DLL_IMPORTEXPORT                    \
    void LamaIO::writeStorageImpl(                      \
            const MatrixStorage<_type>& storage );          \
    \
    template COMMON_DLL_IMPORTEXPORT                    \
    void LamaIO::readStorageImpl(                       \
            MatrixStorage<_type>& storage );                \

SCAI_COMMON_LOOP( SCAI_LAMA_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_LAMA_METHOD_INSTANTIATIONS

}  // lama

}  // scai
