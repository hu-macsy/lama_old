/**
 * @file Communicator.cpp
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
 * @brief Communicator.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Halo.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>
#include <scai/common/Settings.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <locale>
#include <string>
#include <memory>

using namespace std;
using namespace scai::hmemo;
using namespace scai::tasking;

namespace scai
{

namespace dmemo
{

std::ostream& operator<<( std::ostream& stream, const _Communicator::CommunicatorKind& type )
{
    switch ( type )
    {
        case _Communicator::NO :
            stream << "NO";
            break;

        case _Communicator::MPI :
            stream << "MPI";
            break;

        default:
            stream << "CommunicatorKind_" << ( int ) type;
    }

    return stream;
}

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( Communicator::logger, "Communicator" )

CommunicatorPtr Communicator::getCommunicatorPtr( const CommunicatorKind& type )
{
    return create( type );
}

CommunicatorPtr Communicator::getDefaultCommunicatorPtr()
{
    // try MPI communicator for default
    if ( canCreate( MPI ) )
    {
        return create( MPI );
    }

    // if even NO is not availabe an exception is thrown
    return create( NO );
}

CommunicatorPtr Communicator::getCommunicatorPtr()
{
    if ( !CommunicatorStack::empty() )
    {
        return CommunicatorStack::top();
    }

    std::string comm;
    std::locale loc;

    if ( common::Settings::getEnvironment( comm, "SCAI_COMMUNICATOR" ) )
    {
        // comm name not case sensitive, take it upper case
        for ( std::string::iterator p = comm.begin(); comm.end() != p; ++p )
        {
            *p = std::toupper( *p, loc );
        }

        if ( comm == "MPI" )
        {
            return getCommunicatorPtr( MPI );
        }

        if ( comm == "NO" )
        {
            return getCommunicatorPtr( NO );
        }

        COMMON_THROWEXCEPTION( "SCAI_COMMUNICATOR=" << comm << ", unknown communicator type" )
    }

    return getDefaultCommunicatorPtr();
}

const Communicator& Communicator::getCurrent()
{
    return *getCommunicatorPtr();
}

/* -----------------------------------------------------------------------------*/

Communicator::Communicator( const CommunicatorKind& type ) :

    mCommunicatorType( type ),
    mRank( 0 ),
    mSize( 1 ),
    mNodeRank( 0 ),
    mNodeSize( 1 ),
    mSeed( 4711 )
{
    SCAI_LOG_DEBUG( logger, "Communicator constructed, type = " << type )
}

Communicator::~Communicator()
{
    SCAI_LOG_DEBUG( logger, "~Communicator" )
}

bool Communicator::operator==( const Communicator& other ) const
{
    return isEqual( other );
}

bool Communicator::operator!=( const Communicator& other ) const
{
    return !isEqual( other );
}

void Communicator::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "Communicator( type = " << mCommunicatorType << " )";
}

/* -------------------------------------------------------------------------- */

void Communicator::setSizeAndRank( PartitionId size, PartitionId rank )
{
    mSize = size;
    mRank = rank;
    setSeed( mSeed );
}

void Communicator::setSeed( int seed ) const
{
    if ( seed == mSeed )
    {
        return;
    }

    mSeed = seed;

    if ( mSize == 1 )
    {
        // all processors must have the same random numbers
        common::Math::srandom( mSeed );
    }
    else
    {
        // processors must generate different numbers
        common::Math::srandom( mSeed + mRank );
    }
}

/* -------------------------------------------------------------------------- */

void Communicator::factorize2( PartitionId procgrid[2], const double sizeX, const double sizeY ) const
{
    PartitionId usergrid[3];
    getUserProcArray( usergrid );
    // assign partitions to 2d grid so as to minimize contact
    double bestline = 2.0 * ( sizeX + sizeY );
    bool found = false;
    // try all possible factorizations of size
    PartitionId size = getSize();

    for ( PartitionId ipx = 1; ipx <= size; ipx++ )
    {
        if ( size % ipx != 0 )
        {
            continue;
        }

        if ( usergrid[0] && ( usergrid[0] != ipx ) )
        {
            continue;
        }

        PartitionId ipy = size / ipx;

        if ( usergrid[1] && ( usergrid[1] != ipy ) )
        {
            continue;
        }

        double line = sizeX / ipx + sizeY / ipy;

        if ( line < bestline )
        {
            found = true;
            bestline = line;
            procgrid[0] = ipx;
            procgrid[1] = ipy;
            ;
        }
    }

    if ( !found )
    {
        COMMON_THROWEXCEPTION(
            "No processor 2D-grid found for usergrid " << usergrid[0] << " x " << usergrid[1] << ", NP = " << size );
    }

    SCAI_LOG_INFO( logger,
                   "Best processor factorization of size = " << size << ": " << procgrid[0] << " x " << procgrid[1] )
}

/* -------------------------------------------------------------------------- */

void Communicator::factorize3(
    PartitionId procgrid[3],
    const double sizeX,
    const double sizeY,
    const double sizeZ ) const
{
    PartitionId usergrid[3];
    getUserProcArray( usergrid );
    // assign partitions to 3d grid so as to minimize surface area
    double area[3] =
    { sizeX* sizeY, sizeX* sizeZ, sizeY* sizeZ };
    double bestsurf = 2.0 * ( area[0] + area[1] + area[2] );
    // try all possible factorizations of size
    // surface = surface area of a proc sub-domain
    // for 2d, insure ipz = 1
    PartitionId size = getSize();
    bool found = false;

    for ( PartitionId ipx = 1; ipx <= size; ipx++ )
    {
        if ( size % ipx != 0 )
        {
            continue;
        }

        if ( usergrid[0] && ( usergrid[0] != ipx ) )
        {
            continue;
        }

        PartitionId nremain = size / ipx;

        for ( PartitionId ipy = 1; ipy <= nremain; ipy++ )
        {
            if ( nremain % ipy != 0 )
            {
                continue;
            }

            if ( usergrid[1] && ( usergrid[1] != ipy ) )
            {
                continue;
            }

            PartitionId ipz = nremain / ipy;

            if ( usergrid[2] && ( usergrid[2] != ipz ) )
            {
                continue;
            }

            double surf = area[0] / ipx / ipy + area[1] / ipx / ipz + area[2] / ipy / ipz;

            if ( surf < bestsurf )
            {
                found = true;
                bestsurf = surf;
                procgrid[0] = ipx;
                procgrid[1] = ipy;
                procgrid[2] = ipz;
            }
        }
    }

    if ( !found )
    {
        COMMON_THROWEXCEPTION(
            "No processor 3D-grid found for usergrid " << usergrid[0] << " x " << usergrid[1] << " x " << usergrid[2] << ", NP = " << size );
    }

    SCAI_LOG_INFO( logger,
                   "Best processor factorization of size = " << size << ": " << procgrid[0] << " x " << procgrid[1] << " x " << procgrid[2] )
}

/* -------------------------------------------------------------------------- */

void Communicator::getGrid2Rank( PartitionId pos[2], const PartitionId procgrid[2] ) const
{
    SCAI_ASSERT_EQ_ERROR( getSize(), procgrid[0] * procgrid[1], "size mismatch for 2D grid" )
    PartitionId rank = getRank();
    PartitionId size = procgrid[0];
    pos[1] = rank / size;
    rank = rank % size;
    pos[0] = rank;
    SCAI_LOG_INFO( logger,
                   *this << ": is (" << pos[0] << "," << pos[1] << ") of (" << procgrid[0] << "," << procgrid[1] << ")" )
}

/* -------------------------------------------------------------------------- */

void Communicator::getGrid3Rank( PartitionId pos[3], const PartitionId procgrid[3] ) const
{
    SCAI_ASSERT_EQ_ERROR( getSize(), procgrid[0] * procgrid[1] * procgrid[2], "size mismatch for 3D grid" )
    PartitionId rank = getRank();
    PartitionId size = procgrid[0] * procgrid[1];
    pos[2] = rank / size;
    rank = rank % size;
    size = procgrid[0];
    pos[1] = rank / size;
    rank = rank % size;
    pos[0] = rank;
    SCAI_LOG_INFO( logger,
                   *this << ": is (" << pos[0] << "," << pos[1] << "," << pos[2] << ") of (" << procgrid[0] << "," << procgrid[1] << "," << procgrid[2] << ")" )
}

/* -------------------------------------------------------------------------- */

void Communicator::getUserProcArray( PartitionId userProcArray[3] )
{
    std::string npString;

    bool hasNP = common::Settings::getEnvironment( npString, "SCAI_NP" );

    userProcArray[0] = 0;
    userProcArray[1] = 0;
    userProcArray[2] = 0;

    if ( hasNP )
    {
        const std::string delimiters = " x_";

        int offset = 0;
        std::string::size_type lastPos = npString.find_first_not_of( delimiters, 0 );
        // Find first "non-delimiter".
        std::string::size_type pos = npString.find_first_of( delimiters, lastPos );

        while ( std::string::npos != pos || std::string::npos != lastPos )
        {
            // Found a token
            std::istringstream val( npString.substr( lastPos, pos - lastPos ) );

            if ( offset > 2 )
            {
                break; // ignore more than 3 values
            }

            val >> userProcArray[offset++];
            // Skip delimiters.  Note the "not_of"
            lastPos = npString.find_first_not_of( delimiters, pos );
            // Find next "non-delimiter"
            pos = npString.find_first_of( delimiters, lastPos );
        }

        SCAI_LOG_INFO( logger,
                       "SCAI_NP=" << npString << " -> userProcArray " << userProcArray[0] << " x " << userProcArray[1] << " x " << userProcArray[2] )
    }
    else
    {
        SCAI_LOG_INFO( logger, "environment variable SCAI_NP no set" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType Communicator::shift0(
    ValueType targetVals[],
    const IndexType maxTargetSize,
    const ValueType sourceVals[],
    const IndexType sourceSize ) const
{
    SCAI_ASSERT_ERROR( sourceSize <= maxTargetSize, "insufficient size for target array" )

    for ( IndexType i = 0; i < sourceSize; i++ )
    {
        targetVals[i] = sourceVals[i];
    }

    return sourceSize; // effective number of copied values.
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::sumArray( HArray<ValueType>& array ) const
{
    ContextPtr commContext = getCommunicationContext( array );

    SCAI_LOG_INFO( logger, "sumArray<" << common::TypeTraits<ValueType>::id()
                   << " at this context " << *commContext << ", array = " << array )

    SCAI_CONTEXT_ACCESS( commContext );

    IndexType numElems = array.size();

    WriteAccess<ValueType> data( array, commContext );

    // alias of inValues and outValues, sumImpl can deal with it

    sumImpl( data.get(), data.get(), numElems, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::shiftArray(
    HArray<ValueType>& recvArray,
    const HArray<ValueType>& sendArray,
    const int direction ) const
{
    SCAI_ASSERT_UNEQUAL( &recvArray, &sendArray, "send and receive array are same, not allowed for shift" )

    if ( direction % getSize() == 0 )
    {
        // self assignment
        recvArray = sendArray;
        return;
    }

    ContextPtr commContext = getCommunicationContext( sendArray );

    SCAI_LOG_INFO( logger, "shiftArray<" << common::TypeTraits<ValueType>::id()
                   << " at this context " << *commContext << ", sendArray = " << sendArray
                   << ", recvArray = " << recvArray )

    SCAI_CONTEXT_ACCESS( commContext )

    ReadAccess<ValueType> sendData( sendArray, commContext );
    IndexType numSendElems = sendData.size();
    // make recv array large enough to fit for send data
    WriteOnlyAccess<ValueType> recvData( recvArray, commContext, numSendElems );
    // but we are able to receive even more data if array is large enough
    IndexType maxNumRecvElems = recvData.capacity();
    // For shifting of data we use the pure virtual methods implemened by each communicator
    IndexType numRecvElems = shift( recvData.get(), maxNumRecvElems, sendData.get(), numSendElems, direction );
    SCAI_LOG_DEBUG( logger,
                    "shift, direction = " << direction << ", sent " << numSendElems
                    << ", recvd " << numRecvElems << "( max was " << maxNumRecvElems << ")" )
    recvData.resize( numRecvElems ); // take over the size
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::joinArray(
    HArray<ValueType>& globalArray,
    const HArray<ValueType>& localArray ) const
{
    // ToDo: make this routine more efficient, can use gatherAll if sizes are known before

    HArray<ValueType> bufferArray;

    globalArray.clear();

    IndexType currentSize = 0;

    for ( PartitionId p = 0; p < getSize(); ++p )
    {
        if ( p == getRank() )
        {
            // my turn for broadcast

            utilskernel::HArrayUtils::_assign( bufferArray, localArray );
        }

        bcastArray( bufferArray, p );

        // append the local part to the global array

        {
            WriteAccess<ValueType> wGlobal( globalArray );
            ReadAccess<ValueType> rLocal( bufferArray );
            wGlobal.resize( currentSize + rLocal.size() );

            for ( IndexType i = 0; i < rLocal.size(); ++i )
            {
                wGlobal[currentSize + i] = rLocal[i];
            }
        }

        currentSize += bufferArray.size();
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::bcastArray( HArray<ValueType>& array, const IndexType n, const PartitionId root ) const
{
    SCAI_ASSERT_VALID_INDEX_ERROR( root, getSize(), "illegal root for bcast specified" )

    ContextPtr commContext = getCommunicationContext( array );

    SCAI_CONTEXT_ACCESS( commContext )

    if ( root == getRank() )
    {
        SCAI_ASSERT_LE_ERROR( n, array.size(), "bcastArray: root has not enough values" )
        SCAI_LOG_INFO( logger, *this << ": bcast to all other procs, array = " << array )

        ReadAccess<ValueType> rArray( array, commContext );   // WriteAccess might invalidate data

        ValueType* arrayPtr = const_cast<ValueType*>( rArray.get() );

        bcast( arrayPtr, n, root ); // bcast the array
    }
    else
    {
        WriteOnlyAccess<ValueType> wArray( array, commContext, n );
        bcast( wArray.get(), n, root ); // bcast the array
        SCAI_LOG_INFO( logger, *this << ": bcast from root = " << root << ", array = " << array )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::scatterVArray( 
    HArray<ValueType>& recvArray,
    const PartitionId root, 
    const HArray<ValueType>& sendArray, 
    const HArray<IndexType>& sendSizes ) const
{
    // scatter the send sizes so that each processor knows the size of its part to allocate it

    IndexType recvSize;
    ValueType dummy;

    // avoid NULL pointers on non-root processors, caused problems with some MPI implementations

    const IndexType* sizesPtr = &recvSize;
    const ValueType* allValuesPtr = &dummy;

    if ( getRank() == root )
    {
        sizesPtr = hostReadAccess( sendSizes ).get();
        allValuesPtr = hostReadAccess( sendArray ).get();
    }

    scatter( &recvSize, 1, root, sizesPtr );

    auto wMyValues = hostWriteOnlyAccess( recvArray, recvSize );

    scatterV( wMyValues.get(), recvSize, root, allValuesPtr, sizesPtr );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::bcastArray( HArray<ValueType>& array, const PartitionId root ) const
{
    SCAI_ASSERT_VALID_INDEX_ERROR( root, getSize(), "illegal root for bcast specified" )

    ContextPtr commContext = getCommunicationContext( array );

    IndexType n = 0;

    if ( root == getRank() )
    {
        n = array.size();
    }

    bcast( &n, 1, root );

    bcastArray( array, n, root );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* Communicator::shiftAsync(
    HArray<ValueType>& recvArray,
    const HArray<ValueType>& sendArray,
    const int direction ) const
{
    SCAI_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    // ToDo: not quite clear how to deal with asynchronous communication on other devices

    ContextPtr commContext = Context::getHostPtr();

    SCAI_CONTEXT_ACCESS( commContext )

    recvArray.clear(); // do not keep any old data, keep capacities
    WriteAccess<ValueType> recvData( recvArray, commContext );
    ReadAccess<ValueType> sendData( sendArray, commContext );
    IndexType numElems = sendData.size();
    recvData.resize( numElems ); // size should fit at least to keep own data
    // For shifting of data we use the pure virtual methods implemened by each communicator
    // Note: get is the method of the accesses and not of the auto_ptr
    std::unique_ptr<SyncToken> syncToken( shiftAsync( recvData.get(), sendData.get(), numElems, direction ) );
    SCAI_ASSERT_DEBUG( syncToken.get(), "NULL pointer for sync token" )
    // release of accesses are delayed, add routines  in the sync token so they are called at synchonization
    syncToken->pushRoutine( sendData.releaseDelayed() );
    syncToken->pushRoutine( recvData.releaseDelayed() );
    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::updateHalo(
    HArray<ValueType>& haloValues,
    const HArray<ValueType>& localValues,
    const Halo& halo ) const
{
    SCAI_REGION( "Communicator.updateHalo" )
    SCAI_LOG_INFO( logger, *this << ": update halo" )
    const CommunicationPlan& requiredPlan = halo.getRequiredPlan();
    SCAI_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )
    const CommunicationPlan& providesPlan = halo.getProvidesPlan();
    SCAI_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send
    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();
    HArray<ValueType> sendValues( numSendValues ); //!< temporary array for send communication

    utilskernel::HArrayUtils::gather( sendValues, localValues, halo.getProvidesIndexes(), common::BinaryOp::COPY );

    exchangeByPlan( haloValues, requiredPlan, sendValues, providesPlan );
}

/* -------------------------------------------------------------------------- */

static void releaseArray( std::shared_ptr<_HArray> array )
{
    array->clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* Communicator::updateHaloAsync(
    HArray<ValueType>& haloValues,
    const HArray<ValueType>& localValues,
    const Halo& halo ) const
{
    SCAI_REGION( "Communicator.updateHaloAsync" )
    SCAI_LOG_INFO( logger, *this << ": asynchronous update halo" )
    const CommunicationPlan& requiredPlan = halo.getRequiredPlan();
    SCAI_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )
    const CommunicationPlan& providesPlan = halo.getProvidesPlan();
    SCAI_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send
    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();

    std::shared_ptr<HArray<ValueType> > sendValues( new HArray<ValueType>( numSendValues ) );

    // put together the (send) values to provide for other partitions

    utilskernel::HArrayUtils::gather( *sendValues, localValues, halo.getProvidesIndexes(), common::BinaryOp::COPY );

    SyncToken* token( exchangeByPlanAsync( haloValues, requiredPlan, *sendValues, providesPlan ) );

    // Also push the sendValues array to the token so it will be freed after synchronization
    // Note: it is guaranteed that access to sendValues is freed before sendValues
    token->pushRoutine( std::bind( releaseArray, sendValues ) );
    return token;
}

/* -------------------------------------------------------------------------- */
/*  computeOwners                                                             */
/* -------------------------------------------------------------------------- */

void Communicator::computeOwners(
    hmemo::HArray<PartitionId>& owners,
    const Distribution& distribution,
    const hmemo::HArray<IndexType>& requiredIndexes ) const
{
    hmemo::ContextPtr ctx = hmemo::Context::getHostPtr();

    IndexType numIndexes = requiredIndexes.size();

    hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, ctx, numIndexes );
    hmemo::ReadAccess<IndexType> rIndexes( requiredIndexes, ctx );

    computeOwners( wOwners.get(), distribution, rIndexes.get(), numIndexes );
}

/* -------------------------------------------------------------------------- */

void Communicator::computeOwners(
    PartitionId owners[],
    const Distribution& distribution,
    const IndexType requiredIndexes[],
    const IndexType numIndexes ) const
{
    if ( distribution.isReplicated() )
    {
        // there might be multiple owners but we take the only the first one

        const IndexType size = distribution.getGlobalSize();

        for ( IndexType i = 0; i < numIndexes; ++i )
        { 
            SCAI_ASSERT_VALID_INDEX_ERROR( requiredIndexes[i], size, "required index out of range" )
            owners[i] = 0;
        }

        return;
    }

    // Note: this routine is only supported on Host, may change in future releases

    PartitionId rank = getRank();
    PartitionId size = getSize();

    SCAI_LOG_INFO( logger, "need owners for " << numIndexes << " global indexes" )

    SCAI_ASSERT_EQ_ERROR( *this, distribution.getCommunicator(), "illegal communicator for computeOwners" )

    IndexType nonLocal = 0;

    // Check for own ownership. Mark needed Owners. Only exchange requests for unknown indexes.

    for ( IndexType i = 0; i < numIndexes; ++i )
    {
        if ( distribution.isLocal( requiredIndexes[i] ) )
        {
            owners[i] = rank;
        }
        else
        {
            nonLocal++;
            owners[i] = invalidPartition;
        }
    }

    SCAI_LOG_INFO( logger, numIndexes - nonLocal << " Indexes are local. Only need to send " << nonLocal << " values." )
    IndexType receiveSize = max( nonLocal ); // --> pure method call
    SCAI_LOG_DEBUG( logger, "max size of receive buffer is " << receiveSize )
    // Allocate the maximal needed size for the communication buffers
    HArray<IndexType> indexesSendArray( receiveSize );
    HArray<IndexType> indexesReceiveArray( receiveSize );
    HArray<IndexType> ownersSendArray( receiveSize );
    HArray<IndexType> ownersReceiveArray( receiveSize );

    ContextPtr contextPtr = Context::getHostPtr();

    {
        WriteAccess<IndexType> indexesSend( indexesSendArray, contextPtr );
        WriteAccess<IndexType> ownersSend( ownersSendArray, contextPtr );
        nonLocal = 0; // reset, counted again

        for ( IndexType i = 0; i < static_cast<IndexType>( numIndexes ); ++i )
        {
            if ( owners[i] == invalidPartition )
            {
                indexesSend[nonLocal++] = requiredIndexes[i];
            }
        }

        // Important: set owners for full buffer of ownersSend

        for ( IndexType i = 0; i < receiveSize; ++i )
        {
            ownersSend[i] = invalidPartition;
        }
    }
    IndexType ownersSize = invalidIndex;
    IndexType currentSize = nonLocal;
    const int direction = 1; // send to right, recv from left

    for ( PartitionId iProc = 0; iProc < size - 1; ++iProc )
    {
        WriteAccess<IndexType> indexesSend( indexesSendArray, contextPtr );
        WriteAccess<IndexType> indexesReceive( indexesReceiveArray, contextPtr );
        WriteAccess<IndexType> ownersSend( ownersSendArray, contextPtr );
        WriteAccess<IndexType> ownersReceive( ownersReceiveArray, contextPtr );
        SCAI_LOG_DEBUG( logger,
                        *this << " shift: recv " << receiveSize << ", send " << currentSize << ", direction = " << direction )
        // --->   Pure method call
        currentSize = shift( indexesReceive.get(), receiveSize, indexesSend.get(), currentSize, direction );
        SCAI_ASSERT_ERROR( ownersSize == invalidIndex || currentSize == ownersSize, "Communication corrupted." )
        SCAI_LOG_DEBUG( logger, "owners size = " << ownersSize << ", current size = " << currentSize )
        IndexType* indexes = indexesReceive.get();
        IndexType* currentOwners = ownersSend.get();
        SCAI_LOG_DEBUG( logger, "check buffer with " << currentSize << " global indexes whether I am owner" )

        for ( IndexType i = 0; i < currentSize; ++i )
        {
            //TODO there should be a blockwise implementation of isLocal
            SCAI_LOG_TRACE( logger,
                            "check global index " << indexes[i] << " with current owner " << currentOwners[i] << ", is local = " << distribution.isLocal( indexes[i] ) )

            if ( currentOwners[i] == invalidIndex && distribution.isLocal( indexes[i] ) )
            {
                SCAI_LOG_TRACE( logger, *this << ": me is owner of global index " << indexes[i] )
                currentOwners[i] = rank;
            }
        }

        //Release so we can swap the Arrays
        indexesReceive.release();
        indexesSend.release();
        indexesReceiveArray.swap( indexesSendArray );
        SCAI_LOG_DEBUG( logger, *this << ": send array with " << currentSize << " owners to right" )

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            SCAI_LOG_TRACE( logger, *this << " send currentOwner[" << i << "] = " << ownersSend[i] )
        }

        // --->   Pure method call
        ownersSize = shift( ownersReceive.get(), receiveSize, ownersSend.get(), currentSize, direction );
        SCAI_LOG_DEBUG( logger, *this << ": recvd array with " << ownersSize << " owners from left" )

        for ( IndexType i = 0; i < ownersSize; i++ )
        {
            SCAI_LOG_TRACE( logger, *this << ": recv currentOwner[" << i << "] = " << ownersReceive[i] )
        }

        //Release so we can swap the Arrays
        ownersSend.release();
        ownersReceive.release();
        ownersReceiveArray.swap( ownersSendArray );
    }

    WriteAccess<IndexType> ownersSend( ownersSendArray, contextPtr );

    for ( IndexType i = 0; i < nonLocal; ++i )
    {
        SCAI_LOG_TRACE( logger,
                        *this << ": final " << i << " of " << nonLocal << ": " << requiredIndexes[i] << ", owner = " << ownersSend[i] )
    }

    // The Owner Indexes are always passed in the same order, so we can insert them easily.

    IndexType nn = 0;

    for ( IndexType i = 0; i < numIndexes; ++i )
    {
        if ( owners[i] == invalidPartition )
        {
            owners[i] = ownersSend[nn++];
        }
    }

    SCAI_ASSERT_EQUAL( nn, nonLocal, "serious mismatch" )
}

/* -------------------------------------------------------------------------- */

bool Communicator::all( const bool flag ) const
{
    IndexType val = 0; // flag is true

    if ( !flag )
    {
        val = 1;
    }

    IndexType allval = sum( val ); // count flags == false
    SCAI_LOG_DEBUG( logger, "sum( " << val << " ) = " << allval )
    return allval == 0;
}

/* -------------------------------------------------------------------------- */

bool Communicator::any( const bool flag ) const
{
    IndexType val = flag ? 1 : 0; //  1 if flag is true
    IndexType allval = sum( val ); // count flags == false
    SCAI_LOG_DEBUG( logger, "sum( " << val << " ) = " << allval )
    return allval > 0;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Communicator::scan( const ValueType localValue ) const
{
    ValueType scanValue;

    scanImpl( &scanValue, &localValue, 1, common::TypeTraits<ValueType>::stype );

    // scanImpl does an inclusive scan

    return scanValue;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Communicator::scanDefault( ValueType localValue ) const
{
    PartitionId size = getSize();
    PartitionId rank = getRank();
    PartitionId root = 0;

    std::unique_ptr<ValueType[]> allValues( new ValueType[ size ] );

    gather( allValues.get(), 1, root, &localValue );
    bcast ( allValues.get(), size, root );

    ValueType runningSum = 0;

    for ( PartitionId p = 0; p <= rank; ++p )
    {
        runningSum += allValues[p];
    }

    return runningSum;
}

/* -------------------------------------------------------------------------- */

void Communicator::bcast( std::string& val, const PartitionId root ) const
{
    /** Broadcast of a string, done in two steps. */
    bool isRoot = getRank() == root;
    IndexType len = 0;    // IndexType is supported by bcast

    if ( isRoot )
    {
        len = val.length();
    }

    // step 1: broadcast length of string
    bcast( &len, 1, root );
    std::vector<char> buffer( len + 1 );
    char* strptr = &buffer[0];

    if ( isRoot )
    {
        strptr = const_cast<char*>( val.c_str() );
    }

    bcast( strptr, len + 1, root );

    if ( !isRoot )
    {
        val = strptr;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::all2all( ValueType recvValues[], const ValueType sendValues[] ) const
{
    all2allImpl( reinterpret_cast<void*>( recvValues ),
                 reinterpret_cast<const void*>( sendValues ), 
                 common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::all2allv( ValueType* recvBuffer[], const IndexType recvCount[],
                             const ValueType* sendBuffer[], const IndexType sendCount[] ) const
{
    all2allvImpl( reinterpret_cast<void**>( recvBuffer ), recvCount,
                  reinterpret_cast<const void**>( sendBuffer ), sendCount,
                  common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::maxlocDefault( ValueType& val, IndexType& location, const PartitionId root ) const
{
    SCAI_LOG_INFO( logger, *this << ": maxlocDefault<" << common::TypeTraits<ValueType>::id() << ">" )

    ValueType maxVal = max( val );

    IndexType myMaxLocation = invalidIndex;

    if ( maxVal == val )
    {
        myMaxLocation = location;
    }

    std::unique_ptr<IndexType[]> allMaxLocations( new IndexType[ getSize() ] );

    gather( allMaxLocations.get(), 1, root, &myMaxLocation );

    if ( getRank() == root )
    {
        // find first defined location

        for ( PartitionId p = 0; p < getSize(); ++p )
        {
            if ( allMaxLocations[p] != invalidIndex )
            {
                location = allMaxLocations[p];
                SCAI_LOG_DEBUG( logger, *this << ": maxlocDefault location = " << location << " @ " << p )
                break;
            }
        }

        val = maxVal;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::minlocDefault( ValueType& val, IndexType& location, const PartitionId root ) const
{
    SCAI_LOG_INFO( logger, *this << ": minlocDefault<" << common::TypeTraits<ValueType>::id() << ">" )

    ValueType minVal = min( val );

    IndexType myMinLocation = invalidIndex;   // undefined

    if ( minVal == val )
    {
        myMinLocation = location;
    }

    std::unique_ptr<IndexType[]> allMinLocations( new IndexType[ getSize() ] );

    gather( allMinLocations.get(), 1, root, &myMinLocation );

    if ( getRank() == root )
    {
        // find first defined location

        for ( PartitionId p = 0; p < getSize(); ++p )
        {
            if ( allMinLocations[p] != invalidIndex )
            {
                location = allMinLocations[p];
                SCAI_LOG_DEBUG( logger, *this << ": minlocDefault location = " << location << " @ " << p )
                break;
            }
        }

        val = minVal;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::maxloc( ValueType& val, IndexType& location, const PartitionId root ) const
{
    const common::ScalarType vType = common::TypeTraits<ValueType>::stype;
    const common::ScalarType iType = common::TypeTraits<IndexType>::stype;

    if ( supportsLocReduction( vType, iType ) )
    {

        // For the virtual routine maxlocImpl we make sure that val and location are stored contiguously

        struct ValAndLoc
        {
            ValueType val;
            IndexType loc;
        };

        ValAndLoc x;

        x.val = val;
        x.loc = location;

        maxlocImpl( &x.val, &x.loc, root, vType );

        if ( getRank() == root )
        {
            val = x.val;
            location = x.loc;
        }
    }
    else
    {
        maxlocDefault( val, location, root );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::minloc( ValueType& val, IndexType& location, const PartitionId root ) const
{
    const common::ScalarType vType = common::TypeTraits<ValueType>::stype;
    const common::ScalarType iType = common::TypeTraits<IndexType>::stype;

    if ( supportsLocReduction( vType, iType ) )
    {

        // For the virtual routine minlocImpl we make sure that val and location are stored contiguously

        struct ValAndLoc
        {
            ValueType val;
            IndexType loc;
        };

        ValAndLoc x;

        x.val = val;
        x.loc = location;

        minlocImpl( &x.val, &x.loc, root, vType );

        if ( getRank() == root )
        {
            val = x.val;
            location = x.loc;
        }
    }
    else
    {
        minlocDefault( val, location, root );
    }
}

/* -------------------------------------------------------------------------- */

CommunicationPlan Communicator::transpose( const CommunicationPlan& plan ) const
{
    const PartitionId np = getSize();

    // plan might be compressed, so build values again

    std::vector<IndexType> sendSizes( np, 0 );

    for ( PartitionId i = 0; i < plan.size(); ++i )
    {
        const CommunicationPlan::Entry& entry = plan[i];
        sendSizes[entry.partitionId] = entry.quantity;
    }

    std::vector<IndexType> recvSizes( np, 0 );

    // send each processor the number of indexes I require
    // and receive the number of indexes that I have to provide

    all2all( &recvSizes[0], &sendSizes[0] );

    // now we can construct it by quantities

    return CommunicationPlan( recvSizes.data(), np );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::exchangeByPlan(
    hmemo::HArray<ValueType>& recvArray,
    const CommunicationPlan& recvPlan,
    const hmemo::HArray<ValueType>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    SCAI_ASSERT_EQ_ERROR( sendArray.size(), sendPlan.totalQuantity(), "size mismatch" )

    IndexType recvSize = recvPlan.totalQuantity();

    // find a context where data of sendArray can be communicated
    // if possible try to find a context where valid data is available
    // CUDAaware MPI: might give GPU or Host context here

    hmemo::ContextPtr comCtx = getCommunicationContext( sendArray );

    SCAI_LOG_INFO( logger, *this << ": exchangeByPlan<" << common::TypeTraits<ValueType>::id() << ">"
                   << ", send " << sendArray.size() << " values to " << sendPlan.size() << " processors"
                   << ", recv " << recvSize << " values from " << recvPlan.size() << " processors"
                   << ", data at this context " << *comCtx )

    SCAI_CONTEXT_ACCESS( comCtx )

    hmemo::ReadAccess<ValueType> sendData( sendArray, comCtx );
    // Data will be received at the same context where send data is
    hmemo::WriteOnlyAccess<ValueType> recvData( recvArray, comCtx, recvSize );
    exchangeByPlan( recvData.get(), recvPlan, sendData.get(), sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* Communicator::exchangeByPlanAsync(
    hmemo::HArray<ValueType>& recvArray,
    const CommunicationPlan& recvPlan,
    const hmemo::HArray<ValueType>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    SCAI_ASSERT_EQ_ERROR( sendArray.size(), sendPlan.totalQuantity(), "size mismatch" )
    IndexType recvSize = recvPlan.totalQuantity();
    // allocate accesses, SyncToken will take ownership
    hmemo::ContextPtr comCtx = getCommunicationContext( sendArray );
    SCAI_LOG_DEBUG( logger, *this << ": exchangeByPlanAsync, comCtx = " << *comCtx )
    hmemo::ReadAccess<ValueType> sendData( sendArray, comCtx );
    hmemo::WriteOnlyAccess<ValueType> recvData( recvArray, comCtx, recvSize );
    SyncToken* token( exchangeByPlanAsync( recvData.get(), recvPlan, sendData.get(), sendPlan ) );
    // Add the read and write access to the sync token to get it freed after successful wait
    // conversion std::shared_ptr<hmemo::HostWriteAccess<ValueType> > -> std::shared_ptr<BaseAccess> supported
    token->pushRoutine( recvData.releaseDelayed() );
    token->pushRoutine( sendData.releaseDelayed() );
    // return ownership of new created object
    return token;
}

/* -------------------------------------------------------------------------- */

void Communicator::setNodeData()
{
    // Note: here we assume that mSize, mRank are already set correctly

    // routine set mNodeRank and mNodeSize
    // processors with same processor_name are assumed to be on the same node

    int maxNameLength = maxProcessorName();

    std::unique_ptr<char[]> myNodeName( new char[ maxNameLength ] );

    std::unique_ptr<char[]> allNodeNames( new char[ maxNameLength * getSize() ] );

    getProcessorName( myNodeName.get() );

    SCAI_LOG_INFO( logger, "Node name of processor " << mRank << " of " << mSize << ": " << myNodeName.get() )

    memset( allNodeNames.get(), '\0', maxNameLength * mSize * sizeof( char ) );

    // use gather / bcast

    const PartitionId root = 0;

    gather( allNodeNames.get(), maxNameLength, root, myNodeName.get() );
    bcast( allNodeNames.get(), maxNameLength * mSize, root );

    mNodeSize = 0;
    mNodeRank = mSize; // illegal value to verify that it will be set

    const char* ptrAllNodeNames = allNodeNames.get();

    for ( PartitionId i = 0; i < mSize; ++i )
    {
        if ( strcmp( &ptrAllNodeNames[i * maxNameLength], myNodeName.get() ) )
        {
            continue; // processor i is not on same node
        }

        // Processor i is on same node as this processor

        if ( i == mRank )
        {
            mNodeRank = mNodeSize;
        }

        ++mNodeSize;
    }

    SCAI_ASSERT_GT_ERROR( mNodeSize, 0, "Serious problem encountered to get node size" )
    SCAI_ASSERT_LT_ERROR( mNodeRank, mNodeSize, "Serious problem encountered to get node size" )

    SCAI_LOG_INFO( logger, "Processor " << mRank << ": node rank " << mNodeRank << " of " << mNodeSize )
}

/* -------------------------------------------------------------------------- */

#define SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS( _type )             \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    IndexType Communicator::shift0(                                 \
            _type targetVals[],                                     \
            const IndexType maxTargetSize,                          \
            const _type sourceVals[],                               \
            const IndexType sourceSize ) const;                     \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::maxloc(                                      \
            _type& val,                                             \
            IndexType& location,                                    \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::maxlocDefault(                               \
            _type& val,                                             \
            IndexType& location,                                    \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::minloc(                                      \
            _type& val,                                             \
            IndexType& location,                                    \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::minlocDefault(                               \
            _type& val,                                             \
            IndexType& location,                                    \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    _type Communicator::scan( _type val ) const;                    \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    _type Communicator::scanDefault( _type val ) const;             \
                                                                    \
    // instantiate methods for all communicator data types

SCAI_COMMON_LOOP( SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS, SCAI_ALL_TYPES )

#undef SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS

#define SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS( _type )             \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::bcastArray(                                  \
            HArray<_type>& array,                                   \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::bcastArray(                                  \
            HArray<_type>& array,                                   \
            const IndexType n,                                      \
            const PartitionId root ) const;                         \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::shiftArray(                                  \
            HArray<_type>& recvArray,                               \
            const HArray<_type>& sendArray,                         \
            const int direction ) const;                            \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::joinArray(                                   \
            HArray<_type>& globalArray,                             \
            const HArray<_type>& localArray ) const;                \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::scatterVArray(                               \
            HArray<_type>& recvArray,                               \
            const PartitionId root,                                 \
            const HArray<_type>& sendArray,                         \
            const HArray<IndexType>& sendSizes ) const;             \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    SyncToken* Communicator::shiftAsync(                            \
            HArray<_type>& recvArray,                               \
            const HArray<_type>& sendArray,                         \
            const int direction ) const;                            \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::sumArray(                                    \
            HArray<_type>& array ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::updateHalo(                                  \
            HArray<_type>& haloValues,                              \
            const HArray<_type>& localValues,                       \
            const Halo& halo ) const;                               \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::all2all(                                     \
            _type recvValues[],                                     \
            const _type sendValues[] ) const;                       \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::all2allv(                                    \
            _type* recvVal[], const IndexType recvCount[],          \
            const _type* sendVal[],                                 \
            const IndexType sendCount[] ) const;                    \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::exchangeByPlan(                              \
            HArray<_type>& recvArray,                               \
            const CommunicationPlan& recvPlan,                      \
            const HArray<_type>& sendArray,                         \
            const CommunicationPlan& sendPlan ) const;              \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    SyncToken* Communicator::exchangeByPlanAsync(                   \
            HArray<_type>& recvArray,                               \
            const CommunicationPlan& recvPlan,                      \
            const HArray<_type>& sendArray,                         \
            const CommunicationPlan& sendPlan ) const;              \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    SyncToken* Communicator::updateHaloAsync(                       \
            HArray<_type>& haloValues,                              \
            const HArray<_type>& localValues,                       \
            const Halo& halo ) const;

// instantiate communicator methods with Harray only for supported array types

SCAI_COMMON_LOOP( SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_DMEMO_COMMUNICATOR_INSTANTIATIONS

} /* end namespace dmemo */

} /* end namespace scai */
