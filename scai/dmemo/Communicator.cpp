/**
 * @file Communicator.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Communicator.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/Communicator.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Halo.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>
#include <scai/common/Settings.hpp>

using namespace std;
using namespace scai::hmemo;
using namespace scai::tasking;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( Communicator::logger, "Communicator" )

CommunicatorPtr Communicator::getCommunicator( const communicator::CommunicatorKind& type )
{
    SCAI_LOG_TRACE( logger, "Get communicator of type " << type )

    if ( canCreate( type ) )
    {
        return create( type );
    }
    else
    {
        SCAI_LOG_WARN( logger, "could not get communicator " << type << ", take default one" )
        return getCommunicator();
    }
}

CommunicatorPtr Communicator::getCommunicator()
{
    std::string comm;

    if ( common::Settings::getEnvironment( comm, "SCAI_COMMUNICATOR" ) )
    {
        // comm name not case sensitive, take it upper case
        for ( std::string::iterator p = comm.begin(); comm.end() != p; ++p )
        {
            *p = toupper( *p );
        }
        if ( comm == "MPI" ) 
        {
            return getCommunicator( communicator::MPI );
        }
        if ( comm == "GPI" ) 
        {
            return getCommunicator( communicator::GPI );
        }
        if ( comm == "NO" ) 
        {
            return getCommunicator( communicator::NO );
        }
 
        COMMON_THROWEXCEPTION( "SCAI_COMMUNICATOR=" << comm << ", unknown communicator type" )
    }

    // try MPI communicator for default

    if ( canCreate( communicator::MPI ) )
    {   
        return create( communicator::MPI );
    }
    
    // no MPI, try GPI communicator for default

    if ( canCreate( communicator::GPI ) )
    {   
        return create( communicator::GPI );
    }
    
    // if even NO is not availabe an exception is thrown

    return create( communicator::NO );
}

Communicator::Communicator( const communicator::CommunicatorKind& type )
    : mCommunicatorType( type )
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
    stream << "Communicator";
}

void Communicator::factorize2( const double sizeX, const double sizeY, PartitionId procgrid[2] ) const
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

void Communicator::factorize3(
    const double sizeX,
    const double sizeY,
    const double sizeZ,
    PartitionId procgrid[3] ) const
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

void Communicator::getUserProcArray( PartitionId userProcArray[3] )
{
    const char* np4lama = getenv( "LAMA_NP" );

    userProcArray[0] = 0;
    userProcArray[1] = 0;
    userProcArray[2] = 0;

    const std::string delimiters = " x_";

    if ( np4lama )
    {
        std::string str( np4lama );

        int offset = 0;

        std::string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
        // Find first "non-delimiter".
        std::string::size_type pos = str.find_first_of( delimiters, lastPos );

        while ( std::string::npos != pos || std::string::npos != lastPos )
        {
            // Found a token
            std::istringstream val( str.substr( lastPos, pos - lastPos ) );

            if ( offset > 2 )
            {
                break; // ignore more than 3 values
            }

            val >> userProcArray[offset++];

            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of( delimiters, pos );
            // Find next "non-delimiter"
            pos = str.find_first_of( delimiters, lastPos );
        }

        SCAI_LOG_INFO( logger,
                       "LAMA_NP=" << np4lama << " -> userProcArray " << userProcArray[0] << " x " << userProcArray[1] << " x " << userProcArray[2] )
    }
    else
    {
        SCAI_LOG_INFO( logger, "environment variable LAMA_NP no set" )
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
void Communicator::shiftArray(
    HArray<ValueType>& recvArray,
    const HArray<ValueType>& sendArray,
    const int direction ) const
{
    SCAI_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    if ( direction % getSize() == 0 )
    {
        // self assignment

        recvArray = sendArray;
        return;
    }

    ContextPtr commContext = getCommunicationContext( sendArray );

    SCAI_LOG_DEBUG( logger,
                    "shiftArray at this context " << *commContext << ", sendArray = " << sendArray << ", recvArray = " << recvArray )

    ReadAccess<ValueType> sendData( sendArray, commContext );

    IndexType numSendElems = sendData.size();

    // make recv array large enough to fit for send data

    WriteOnlyAccess<ValueType> recvData( recvArray, commContext, numSendElems );

    // but we are able to receive even more data if array is large enough

    IndexType maxNumRecvElems = recvData.capacity();

    // For shifting of data we use the pure virtual methods implemened by each communicator

    IndexType numRecvElems = shiftData( recvData.get(), maxNumRecvElems, sendData.get(), numSendElems, direction );

    SCAI_LOG_INFO( logger,
                   "shift, direction = " << direction << ", sent " << numSendElems << ", recvd " << numRecvElems << "( max was " << maxNumRecvElems << ")" )

    recvData.resize( numRecvElems ); // take over the size
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* Communicator::shiftAsync(
    HArray<ValueType>& recvArray,
    const HArray<ValueType>& sendArray,
    const int direction ) const
{
    SCAI_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

    recvArray.clear(); // do not keep any old data, keep capacities

    WriteAccess<ValueType> recvData( recvArray, contextPtr );
    ReadAccess<ValueType> sendData( sendArray, contextPtr );

    IndexType numElems = sendData.size();

    recvData.resize( numElems ); // size should fit at least to keep own data

    // For shifting of data we use the pure virtual methods implemened by each communicator
    // Note: get is the method of the accesses and not of the auto_ptr

    common::unique_ptr<SyncToken> syncToken( shiftDataAsync( recvData.get(), sendData.get(), numElems, direction ) );

    SCAI_ASSERT_DEBUG( syncToken.get(), "NULL pointer for sync token" )

    // release of accesses are delayed, add routines  in the sync token so they are called at synchonization

    syncToken->pushRoutine( sendData.releaseDelayed() );
    syncToken->pushRoutine( recvData.releaseDelayed() );

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::updateHalo(
    HArray<ValueType> &haloValues,
    const HArray<ValueType>& localValues,
    const Halo& halo ) const
{
    SCAI_REGION( "Communicator.updateHalo" )

    SCAI_LOG_INFO( logger, *this << ": update halo" )

    const CommunicationPlan& requiredPlan = halo.getRequiredPlan();

    SCAI_ASSERT_ERROR( requiredPlan.allocated(), "Required plan in Halo not allocated" )
    SCAI_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )

    const CommunicationPlan& providesPlan = halo.getProvidesPlan();

    SCAI_ASSERT_ERROR( providesPlan.allocated(), "Provides plan in Halo not allocated" )
    SCAI_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send

    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();

    HArray<ValueType> sendValues( numSendValues ); //!< temporary array for send communication

    // TODO: HArrayUtils::gather( sendValues, localValues, halo.getProvidesIndexes() );

    {
        ReadAccess<ValueType> local( localValues );
        ReadAccess<IndexType> ind ( halo.getProvidesIndexes() );

        IndexType nValues = ind.size();

        WriteOnlyAccess<ValueType> send ( sendValues, nValues );

        for ( IndexType i = 0; i < nValues; ++i )
        {
            send[i] = local[ ind[i] ];
        }
    }

    exchangeByPlan( haloValues, requiredPlan, sendValues, providesPlan );
}

/* -------------------------------------------------------------------------- */

static void releaseArray( common::shared_ptr<_HArray> array )
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

    SCAI_ASSERT_ERROR( requiredPlan.allocated(), "Required plan in Halo not allocated" )
    SCAI_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )

    const CommunicationPlan& providesPlan = halo.getProvidesPlan();

    SCAI_ASSERT_ERROR( providesPlan.allocated(), "Provides plan in Halo not allocated" )
    SCAI_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send

    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();

    common::shared_ptr<HArray<ValueType> > sendValues( new HArray<ValueType>( numSendValues ) );

    // put together the (send) values to provide for other partitions

    {
        ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

        WriteAccess<ValueType> sendData( *sendValues, contextPtr );
        ReadAccess<ValueType> localData( localValues, contextPtr );
        ReadAccess<IndexType> sendIndexes( halo.getProvidesIndexes(), contextPtr );

        // currently supported only on host

        for ( IndexType i = 0; i < numSendValues; i++ )
        {
            sendData[i] = localData[sendIndexes[i]];
        }
    }

    SyncToken* token( exchangeByPlanAsync( haloValues, requiredPlan, *sendValues, providesPlan ) );

    // Also push the sendValues array to the token so it will be freed after synchronization

    // Note: it is guaranteed that access to sendValues is freed before sendValues

    token->pushRoutine( common::bind( releaseArray, sendValues ) );

    return token;
}

/* -------------------------------------------------------------------------- */
/*  computeOwners                                                             */
/* -------------------------------------------------------------------------- */

void Communicator::computeOwners(
    PartitionId owners[],
    const Distribution& distribution,
    const IndexType requiredIndexes[],
    const IndexType nIndexes ) const
{
    PartitionId rank = getRank();
    PartitionId size = getSize();

    SCAI_LOG_DEBUG( logger, "need owners for " << nIndexes << " global indexes" )

    if ( distribution.getCommunicator() != *this )
    {
        COMMON_THROWEXCEPTION( "The distribution has a different Communicator." )
    }

    IndexType nonLocal = 0;

    // Check for own ownership. Mark needed Owners. Only exchange requests for unknown indexes.
    for ( IndexType i = 0; i < nIndexes; ++i )
    {
        if ( distribution.isLocal( requiredIndexes[i] ) )
        {
            owners[i] = rank;
        }
        else
        {
            nonLocal++;
            owners[i] = nIndex;
        }
    }

    SCAI_LOG_DEBUG( logger, nIndexes - nonLocal << " Indexes are local. Only need to send " << nonLocal << " values." )
    IndexType receiveSize = max( nonLocal ); // --> pure method call
    SCAI_LOG_DEBUG( logger, "max size of receive buffer is " << receiveSize )

    // Allocate the maximal needed size for the communication buffers

    HArray<IndexType> indexesSendArray( receiveSize );
    HArray<IndexType> indexesReceiveArray( receiveSize );
    HArray<IndexType> ownersSendArray( receiveSize );
    HArray<IndexType> ownersReceiveArray( receiveSize );

    ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

    {
        WriteAccess<IndexType> indexesSend( indexesSendArray, contextPtr );
        WriteAccess<IndexType> ownersSend( ownersSendArray, contextPtr );

        nonLocal = 0; // reset, counted again

        for ( IndexType i = 0; i < static_cast<IndexType>( nIndexes ); ++i )
        {
            if ( owners[i] == nPartition )
            {
                indexesSend[nonLocal++] = requiredIndexes[i];
            }
        }

        // Important: set owners for full buffer of ownersSend

        for ( IndexType i = 0; i < receiveSize; ++i )
        {
            ownersSend[i] = nIndex;
        }
    }

    IndexType ownersSize = nIndex;
    IndexType currentSize = nonLocal;

    const int direction = 1; // send to right, recv from left

    for ( int iProc = 0; iProc < size - 1; ++iProc )
    {
        WriteAccess<IndexType> indexesSend( indexesSendArray, contextPtr );
        WriteAccess<IndexType> indexesReceive( indexesReceiveArray, contextPtr );
        WriteAccess<IndexType> ownersSend( ownersSendArray, contextPtr );
        WriteAccess<IndexType> ownersReceive( ownersReceiveArray, contextPtr );
        SCAI_LOG_DEBUG( logger,
                        *this << " shift: recv " << receiveSize << ", send " << currentSize << ", direction = " << direction )

        // --->   Pure method call

        currentSize = shiftData( indexesReceive.get(), receiveSize, indexesSend.get(), currentSize, direction );

        SCAI_ASSERT_ERROR( ownersSize == nIndex || currentSize == ownersSize, "Communication corrupted." )

        SCAI_LOG_DEBUG( logger, "owners size = " << ownersSize << ", current size = " << currentSize )
        IndexType* indexes = indexesReceive.get();
        IndexType* currentOwners = ownersSend.get();
        SCAI_LOG_DEBUG( logger, "check buffer with " << currentSize << " global indexes whether I am owner" )

        for ( IndexType i = 0; i < currentSize; ++i )
        {
            //TODO there should be a blockwise implementation of isLocal
            SCAI_LOG_TRACE( logger,
                            "check global index " << indexes[i] << " with current owner " << currentOwners[i] << ", is local = " << distribution.isLocal( indexes[i] ) )

            if ( currentOwners[i] == nIndex && distribution.isLocal( indexes[i] ) )
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
        ownersSize = shiftData( ownersReceive.get(), receiveSize, ownersSend.get(), currentSize, direction );

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

    for ( int i = 0; i < nonLocal; ++i )
    {
        SCAI_LOG_TRACE( logger,
                        *this << ": final " << i << " of " << nonLocal << ": " << requiredIndexes[i] << ", owner = " << ownersSend[i] )
    }

    // The Owner Indexes are always passed in the same order, so we can insert them easily.

    int nn = 0;

    for ( IndexType i = 0; i < nIndexes; ++i )
    {
        if ( owners[i] == nIndex )
        {
            owners[i] = ownersSend[nn++];

            //TODO is this usefull for the speed ?
            if ( nn == nonLocal )
            {
                break;
            }
        }
    }
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

// Instantiation of template methods for the supported types

#define LAMA_COMMUNICATOR_INSTANTIATIONS( z, I, _ )                 \
    \
    template COMMON_DLL_IMPORTEXPORT                                \
    IndexType Communicator::shift0(                                 \
            ARRAY_TYPE##I targetVals[],                             \
            const IndexType maxTargetSize,                          \
            const ARRAY_TYPE##I sourceVals[],                       \
            const IndexType sourceSize ) const;                     \
    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::shiftArray(                                  \
            HArray<ARRAY_TYPE##I>& recvArray,                    \
            const HArray<ARRAY_TYPE##I>& sendArray,              \
            const int direction ) const;                            \
    \
    template COMMON_DLL_IMPORTEXPORT                                \
    SyncToken* Communicator::shiftAsync(                            \
            HArray<ARRAY_TYPE##I>& recvArray,                    \
            const HArray<ARRAY_TYPE##I>& sendArray,              \
            const int direction ) const;                            \
    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void Communicator::updateHalo(                                  \
            HArray<ARRAY_TYPE##I>& haloValues,                   \
            const HArray<ARRAY_TYPE##I>& localValues,            \
            const Halo& halo ) const;                               \
    \
    template COMMON_DLL_IMPORTEXPORT                                \
    SyncToken* Communicator::updateHaloAsync(                       \
            HArray<ARRAY_TYPE##I>& haloValues,                   \
            const HArray<ARRAY_TYPE##I>& localValues,            \
            const Halo& halo ) const;                               \
     

// instantiate methods for all supported data types

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_COMMUNICATOR_INSTANTIATIONS, _ )

#undef LAMA_COMMUNICATOR_INSTANTIATIONS

} /* end namespace dmemo */

} /* end namespace scai */
