/**
 * @file _Matrix.cpp
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
 * @brief Implementation of methods for all kind of matrices.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/matrix/_Matrix.hpp>

// local library
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/io/PartitionIO.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Settings.hpp>

namespace scai
{

using namespace common;
using namespace dmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( _Matrix::logger, "Matrix" )

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a matrix                                                          */
/* ---------------------------------------------------------------------------------------*/

_Matrix* _Matrix::getMatrix( const Format format, const common::ScalarType valueType )
{
    MatrixCreateKeyType mattype( format, valueType );
    return _Matrix::create( mattype );
}

/* ----------------------------------------------------------------------- */

_Matrix::_Matrix( const _Matrix& other ) :

    Distributed( other ),
    mColDistribution( other.mColDistribution ),
    mCommunicationKind( other.mCommunicationKind )
{
    SCAI_LOG_INFO( logger, "Creating copy of " << other << " with same distributions." )
}

/* ----------------------------------------------------------------------- */

_Matrix::_Matrix( const _Matrix& other, DistributionPtr rowDist, DistributionPtr colDist ) :

    Distributed( rowDist ),
    mColDistribution( colDist ),
    mCommunicationKind( other.mCommunicationKind )
{
    // Very important: here we check that new distributions fit the matrix
    checkSettings();
    SCAI_LOG_INFO( logger,
                   "Creating copy of " << other << " with new distributions: " << "row = " << getDistribution() << ", col = " << getColDistribution() )
}

/* ----------------------------------------------------------------------- */

_Matrix::_Matrix( const IndexType numRows, const IndexType numColumns ) :

    Distributed( DistributionPtr( new NoDistribution( numRows ) ) ),
    mColDistribution( DistributionPtr( new NoDistribution( numColumns ) ) )
{
    setDefaultKind();
    SCAI_LOG_INFO( logger, "Creating a replicated _Matrix of size " << numRows << " x " << numColumns )
}

/* ----------------------------------------------------------------------- */

void _Matrix::setIdentity( const IndexType n )
{
    // take replicated distribution and use pure method
    setIdentity( DistributionPtr( new NoDistribution( n ) ) );
}

/* ----------------------------------------------------------------------- */

void _Matrix::checkSettings() const
{
    if ( !mColDistribution )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for column distribution" )
    }
}

/* ----------------------------------------------------------------------- */

_Matrix::_Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
    : Distributed( rowDistribution )
{
    setDistributedMatrix( rowDistribution, colDistribution );
    setDefaultKind();
    SCAI_LOG_INFO( logger,
                   "Construct a _Matrix of size " << getNumRows() << " x " << getNumColumns() 
                    << " with the distribution " << getDistribution() )
}

_Matrix::_Matrix( DistributionPtr distribution )
    : Distributed( distribution )
{
    setDistributedMatrix( distribution, distribution );
    SCAI_LOG_INFO( logger,
                   "Construct a square _Matrix of size " << getNumRows() << " x " << getNumColumns() 
                   << " with the row/col distribution " << getDistribution() )
}

_Matrix::_Matrix() : 

    Distributed( DistributionPtr( new NoDistribution( 0 ) ) ), 
    mColDistribution( DistributionPtr( new NoDistribution( 0 ) ) )
{
    setDefaultKind();
}

_Matrix::~_Matrix()
{
    SCAI_LOG_DEBUG( logger, "~Matrix" )
}

SyncKind _Matrix::getDefaultSyncKind()
{
    static bool computed = false;

    static SyncKind syncKind = SyncKind::ASYNC_COMM;

    if ( !computed )
    {
        int kind = 1;

        common::Settings::getEnvironment( kind, "SCAI_ASYNCHRONOUS" );

        if ( kind == 0 )
        {
            syncKind = SyncKind::SYNCHRONOUS;
        }
        else if ( kind == 2 )
        {
            syncKind = SyncKind::ASYNC_LOCAL;
        }
    }

    return syncKind;
}

void _Matrix::setDefaultKind()
{
    mCommunicationKind = getDefaultSyncKind();
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::setDistributedMatrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_ASSERT_ERROR( rowDistribution, "NULL row distribution for matrix not allowed" )
    SCAI_ASSERT_ERROR( colDistribution, "NULL column distribution for matrix not allowed" )
    setDistributionPtr( rowDistribution );
    mColDistribution = colDistribution;
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::moveImpl( _Matrix&& other )
{
    auto nullSize = std::make_shared<NoDistribution>( 0 );

    setDistributionPtr( other.getRowDistributionPtr() );
    mColDistribution = other.mColDistribution;

    other.setDistributionPtr( nullSize );
    other.mColDistribution = nullSize;
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::setReplicatedMatrix( const IndexType numRows, const IndexType numColumns )
{
    DistributionPtr rowDist( new NoDistribution( numRows ) );

    if ( numRows == numColumns )
    {
        setDistributedMatrix( rowDist, rowDist );
    }
    else
    {
        setDistributedMatrix( rowDist, DistributionPtr( new NoDistribution( numColumns ) ) );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::setCommunicationKind( SyncKind communicationKind )
{
    mCommunicationKind = communicationKind;
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::inheritAttributes( const _Matrix& other )
{
    setCommunicationKind( other.getCommunicationKind() );
    setContextPtr( other.getContextPtr() );
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeAt( std::ostream& stream ) const
{
    stream << "Matrix(" << getNumRows() << "x" << getNumColumns() << ")";
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::swapMatrix( _Matrix& other )
{
    Distributed::swap( other );
    std::swap( mColDistribution, other.mColDistribution );
}

/* ---------------------------------------------------------------------------------*/

double _Matrix::getSparsityRate() const
{
    return ( double ) getNumValues() / getNumRows() / getNumColumns();
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToSingleFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType /* = UNKNOWN for DEFAULT */,
    const common::ScalarType indexType /* = UNKNOWN for DEFAULT */,
    const FileIO::FileMode fileMode /* = DEFAULT_MODE */ ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    if ( getDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        if ( comm->getRank() == 0 )
        {
            getLocalStorage().writeToFile( fileName, fileType, dataType, indexType, fileMode );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file
        comm->synchronize();
    }
    else
    {
        DistributionPtr rowDist( new NoDistribution( getNumRows() ) );
        DistributionPtr colDist( new NoDistribution( getNumColumns() ) );
        std::unique_ptr<_Matrix> repM( copy( rowDist, colDist ) );
        repM->writeToSingleFile( fileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToPartitionedFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType /* = UNKNOWN for DEFAULT */,
    const common::ScalarType indexType /* = UNKNOWN for DEFAULT */,
    const FileIO::FileMode fileMode /* = DEFAULT_MODE */ ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    if ( getColDistribution().isReplicated() )
    {
        // each processor writes its partition to a file with unique name

        getLocalStorage().writeToFile( fileName, fileType, dataType, indexType, fileMode );
    }
    else
    {
        DistributionPtr colDist( new NoDistribution( getNumColumns() ) );
        std::unique_ptr<_Matrix> repM( copy( getRowDistributionPtr(), colDist ) );
        repM->writeToPartitionedFile( fileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType,
    const common::ScalarType indexType,
    const FileIO::FileMode fileMode ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    std::string newFileName = fileName;

    bool isPartitioned;

    const Communicator& comm = getRowDistribution().getCommunicator();

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, comm );

    if ( !isPartitioned )
    {
        writeToSingleFile( newFileName, fileType, dataType, indexType, fileMode );
    }
    else
    {
        writeToPartitionedFile( newFileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromSingleFile( const std::string& fileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const PartitionId MASTER = 0;
    const PartitionId myRank = comm->getRank();

    // this is a bit tricky stuff, but it avoids an additional copy from storage -> matrix

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    IndexType dims[2];

    if ( myRank == MASTER )
    {
        localMatrix.readFromFile( fileName );

        dims[0] = localMatrix.getNumRows();
        dims[1] = localMatrix.getNumColumns();
    }

    comm->bcast( dims, 2, MASTER );

    if ( myRank != MASTER )
    {
        IndexType localNumRows = 0;
        localMatrix.allocate( localNumRows, dims[1] );
    }

    DistributionPtr rowDist( new SingleDistribution( dims[0], comm, MASTER ) );
    DistributionPtr colDist( new NoDistribution( dims[1] ) );

    // works fine as assign can deal with alias, i.e. localMatrix und getLocalStorage() are same

    SCAI_LOG_DEBUG( logger, *comm << ": assign local storage " << localMatrix );

    assignLocal( localMatrix, rowDist );
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromSingleFile( const std::string& fileName, const DistributionPtr distribution )
{
    if ( distribution.get() == NULL )
    {
        readFromSingleFile( fileName );
        return;
    }

    // dist must be block distributed, not checked again here

    const IndexType n = distribution->getBlockDistributionSize();

    if ( n == invalidIndex )
    {
        readFromSingleFile( fileName );
        redistribute( distribution, getColDistributionPtr() );
        return;
    }

    IndexType first = 0;

    if ( n > 0 )
    {
        first = distribution->local2global( 0 );   // first global index
    }

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    bool error = false;

    try
    {
        localMatrix.readFromFile( fileName, first, n );
    }
    catch ( Exception& ex )
    {
        SCAI_LOG_ERROR( logger, ex.what() )
        error = true;
    }

    error = distribution->getCommunicator().any( error );

    if ( error )
    {
        COMMON_THROWEXCEPTION( "readFromSingleFile failed." )
    }

    // ToDo: what happens if num columns is not same

    assignLocal( localMatrix, distribution );
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromPartitionedFile( const std::string& myPartitionFileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // this is a bit tricky stuff, but it avoids an additional copy from storage -> matrix

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    bool errorFlag = false;

    IndexType localSize = 0;

    try
    {
        localMatrix.readFromFile( myPartitionFileName );

        localSize = localMatrix.getNumRows();

    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_ERROR( logger, *comm << ": failed to read " << myPartitionFileName << ": " << e.what() )
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "error reading partitioned matrix" )
    }

    // We assume a general block distribution

    IndexType globalSize = comm->sum( localSize );

    DistributionPtr rowDist( new GenBlockDistribution( globalSize, localSize, comm ) );

    // make sure that all processors have the same number of columns

    IndexType numColumns = comm->max( localMatrix.getNumColumns() );

    // for consistency we have to reset the number of columns in each stroage

    localMatrix.resetNumColumns( numColumns );

    assignLocal( localMatrix, rowDist );
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::resetRowDistribution( DistributionPtr newDist )
{
    SCAI_ASSERT_EQ_ERROR( getNumRows(), newDist->getGlobalSize(), "global size mismatch" )

    const _MatrixStorage& localMatrix = getLocalStorage();

    SCAI_ASSERT_EQ_ERROR( localMatrix.getNumRows(), newDist->getLocalSize(), "local size mismatch" );

    setDistributionPtr( newDist );
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromFile( const std::string& matrixFileName, const std::string& distributionFileName )
{
    if ( distributionFileName == "BLOCK" )
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        DistributionPtr rowDist;

        // for a single file we set a BlockDistribution

        if ( matrixFileName.find( "%r" ) == std::string::npos )
        {
            PartitionId root = 0;

            IndexType numRows = invalidIndex;

            if ( comm->getRank() == root )
            {
                numRows = FileIO::getStorageSize( matrixFileName );
            }

            comm->bcast( &numRows, 1, root );

            rowDist.reset( new BlockDistribution( numRows, comm ) );
        }

        readFromFile( matrixFileName, rowDist );
    }
    else
    {
        // read the distribution

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        DistributionPtr rowDist = PartitionIO::readDistribution( distributionFileName, comm );

        readFromFile( matrixFileName, rowDist );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromFile( const std::string& fileName, DistributionPtr rowDist )
{
    SCAI_LOG_INFO( logger,
                   *this << ": readFromFile( " << fileName << " )" )

    std::string newFileName = fileName;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // take default

    if ( rowDist.get() )
    {
        comm = rowDist->getCommunicatorPtr();
    }

    bool isPartitioned;

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, *comm );

    SCAI_LOG_INFO( logger, *comm << ": _Matrix.readFromFile ( " << fileName << " ) -> read "
                   << newFileName << ", partitioned = " << isPartitioned );

    if ( !isPartitioned )
    {
        readFromSingleFile( newFileName, rowDist );
    }
    else
    {
        readFromPartitionedFile( newFileName );

        if ( rowDist.get() )
        {
            resetRowDistribution( rowDist );
        }
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::redistribute( const dmemo::Redistributor& redistributor )
{
    if ( getColDistribution().isReplicated() ) 
    {
        redistribute( redistributor, getColDistributionPtr() );
    }
    else if ( getColDistribution() == getRowDistribution() )
    {
        redistribute( redistributor, redistributor.getTargetDistributionPtr() );
    }
    else
    {
        COMMON_THROWEXCEPTION( "redistribute: no new column distribution" )
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::checkLocalStorageSizes( const _MatrixStorage& localStorage, const Distribution& rowDist )
{
    // make some 'global' checks to verify correct sizes on all processors

    const Communicator& comm = rowDist.getCommunicator();

    IndexType maxColumns = comm.max( localStorage.getNumColumns() );

    bool okay = true;

    if ( localStorage.getNumRows() != rowDist.getLocalSize() )
    {
        SCAI_LOG_ERROR( logger, comm << ": #rows of local storage " << localStorage.getNumRows()
                                << " does not match local size of row dist = " << rowDist )
        okay = false;
    }

    if ( localStorage.getNumColumns() != maxColumns )
    {
        SCAI_LOG_ERROR( logger, comm << ": #columns of local storage " << localStorage.getNumColumns()
                               << " must be same on all processors ( max = " << maxColumns << " )" )
        okay = false;
    }

    okay = comm.all( okay );

    if ( !okay )
    { 
        COMMON_THROWEXCEPTION( "Constructor of sparse matrix by local storages failed" )
    } 
}

/* ---------------------------------------------------------------------------------*/

_Matrix* _Matrix::copy( DistributionPtr rowDistribution, DistributionPtr colDistribution ) const
{
    // simple default implementation that works for each matrix
    std::unique_ptr<_Matrix> rep( copy() );
    // unique_ptr guarantees that data is freed if redistribute fails for any reason
    rep->redistribute( rowDistribution, colDistribution );
    return rep.release();
}

MatrixCreateKeyType _Matrix::getCreateValue() const
{
    return MatrixCreateKeyType( getFormat(), getValueType() );
}

} /* end namespace lama */

} /* end namespace scai */
