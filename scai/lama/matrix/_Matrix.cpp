/**
 * @file _Matrix.cpp
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
#include <scai/dmemo/RedistributePlan.hpp>

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

void _Matrix::writeToFile(
    const std::string& fileName,
    const FileMode fileMode,
    const common::ScalarType dataType,
    const common::ScalarType indexType ) const
{
    SCAI_LOG_INFO( logger, *this << ": writeToFile " << fileName )

    std::string suffix = FileIO::getSuffix( fileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> file( FileIO::create( suffix ) );

        file->setMode( fileMode );
        file->setDataType( dataType );
        file->setIndexType( indexType );

        bool appendMode = false;
        common::Settings::getEnvironment( appendMode, "SCAI_IO_APPEND" );
        file->open( fileName.c_str(), appendMode ? "a" : "w" );

        SCAI_LOG_DEBUG( logger, *this << ": writeToFile uses opened " << *file )
        writeToFile( *file );
        file->close();
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToFileBlocked( FileIO& file ) const
{
    // this matrix must be block distributed onto the processors assigned to the file

    CommunicatorPtr comm = file.getCommunicatorPtr();

    DistributionPtr dist = getDistributionPtr();

    SCAI_LOG_DEBUG( logger, "matrix distribution = " << *dist << ", file comm = " << *comm )

    bool isBlockDistributed = dist->isBlockDistributed( comm );

    if ( isBlockDistributed && getColDistribution().isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": independent/collective write, write local storage: " << getLocalStorage() )
        file.writeStorage( getLocalStorage() );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << ": independent/collective write, redistribution to block dist required" )

        auto rowDist = isBlockDistributed ? dist : dist->toBlockDistribution( comm );
        auto colDist = noDistribution( getNumColumns() );
        std::unique_ptr<_Matrix> matBlockDistributed( copyRedistributed( rowDist, colDist ) );
        matBlockDistributed->writeToFile( file );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToFileSingle( FileIO& file ) const
{
    // MASTER mode: only master process 0 writes to file the whole storage

    DistributionPtr dist = getDistributionPtr();

    CommunicatorPtr comm = file.getCommunicatorPtr();

    bool isMasterDistributed = dist->isMasterDistributed( comm );

    if ( isMasterDistributed && getColDistribution().isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": single mode, only MASTER writes it" )

        if ( comm->getRank() == 0 )
        {
            file.writeStorage( getLocalStorage() );
        }
 
        // add a barrier ??, might be better to have it on the close that is called by all processors
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << ": single mode, replicate it" )

        auto rowDist = isMasterDistributed ? dist : dist->toMasterDistribution( comm );
        auto colDist = noDistribution( getNumColumns() );
        std::unique_ptr<_Matrix> matSingle( copyRedistributed( rowDist, colDist ) );
        matSingle->writeToFile( file );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::writeToFile( FileIO& file ) const
{
    SCAI_LOG_INFO( logger, *this << ": writeToFile( file = " << file << " )" )

    if ( file.getDistributedIOMode() == DistributedIOMode::MASTER )
    {
        // master processor will write the whole data

        writeToFileSingle( file );
    }
    else
    {
        // all processors will write the 'block-distributed' data

        writeToFileBlocked( file );
    }
}

/* ---------------------------------------------------------------------------------*/

void _Matrix::readFromFile( const std::string& fileName )
{
    std::string suffix = FileIO::getSuffix( fileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> file( FileIO::create( suffix ) );

        file->open( fileName.c_str(), "r" );
        readFromFile( *file );
        file->close();
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Matrix::readFromFile( const std::string& fileName, DistributionPtr distribution )
{
    // ToDo: distribution might be used for collective mode if it is a block distribution
    readFromFile( fileName );
    redistribute( distribution, getColDistributionPtr() );
}

/* ---------------------------------------------------------------------------------------*/

void _Matrix::readFromFile( FileIO& file )
{
    // Note: in contrary to _Vector we can write it here for all matrix classes by using getLocalStorage()

    SCAI_LOG_INFO( logger, "read matrix from file: " << file )

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    auto comm = dmemo::Communicator::getCommunicatorPtr();   // current communicator

    if ( file.getDistributedIOMode() == DistributedIOMode::MASTER )
    {
        const PartitionId MASTER = 0;
        const PartitionId myRank = comm->getRank();

        IndexType globalDims[2];  // [ #rows, #cols ] 

        if ( myRank == MASTER )
        {
            SCAI_LOG_INFO( logger, *comm << ": read local array from file " << file )
            file.readStorage( localMatrix );
            globalDims[0] = localMatrix.getNumRows(); 
            globalDims[1] = localMatrix.getNumColumns(); 
        }

        comm->bcast( globalDims, 2, MASTER );

        if ( myRank != MASTER )
        {
            IndexType localNumRows = 0;
            localMatrix.allocate( localNumRows, globalDims[1] );
        }

        assignLocal( localMatrix, singleDistribution( globalDims[0], comm, MASTER ) );
    }
    else
    {
        // INDEPENDENT or COLLECTIVE 

        file.readStorage( localMatrix );

        // make sure that all processors have the same number of columns

        IndexType numColumns = comm->max( localMatrix.getNumColumns() );

        // for consistency we have to reset the number of columns in each stroage

        localMatrix.resetNumColumns( numColumns );

        assignLocal( localMatrix, genBlockDistributionBySize( localMatrix.getNumRows(), comm ) );
    }
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

void _Matrix::redistribute( const dmemo::RedistributePlan& redistributor )
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

void _Matrix::replicate()
{
    redistribute( dmemo::noDistribution( getNumRows() ), dmemo::noDistribution( getNumColumns() ) );
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

_Matrix* _Matrix::copyRedistributed( DistributionPtr rowDistribution, DistributionPtr colDistribution ) const
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
