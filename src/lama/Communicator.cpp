/**
 * @file Communicator.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @author Jiri Kraus
 * @date 23.02.2011
 * $Id$
 */

// hpp
#include <lama/Communicator.hpp>

// others
#include <lama/LAMAArrayUtils.hpp>
#include <lama/NoSyncToken.hpp>

#include <lama/distribution/Distribution.hpp>
#include <lama/distribution/Halo.hpp>

// tracing
#include <lama/tracing.hpp>

using namespace std;

namespace lama
{

LAMA_LOG_DEF_LOGGER( Communicator::logger, "Communicator" )

Communicator::Communicator( const std::string& type )
    : mCommunicatorType( type )
{
    LAMA_LOG_DEBUG( logger, "Communicator constructed, type = " << type )
}

Communicator::~Communicator()
{
    LAMA_LOG_DEBUG( logger, "~Communicator" )
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
        LAMA_THROWEXCEPTION(
            "No processor 2D-grid found for usergrid " << usergrid[0] << " x " << usergrid[1] << ", NP = " << size );
    }

    LAMA_LOG_INFO( logger,
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
    { sizeX * sizeY, sizeX * sizeZ, sizeY * sizeZ };

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
        LAMA_THROWEXCEPTION(
            "No processor 3D-grid found for usergrid " << usergrid[0] << " x " << usergrid[1] << " x " << usergrid[2] << ", NP = " << size );
    }

    LAMA_LOG_INFO( logger,
                   "Best processor factorization of size = " << size << ": " << procgrid[0] << " x " << procgrid[1] << " x " << procgrid[2] )
}

void Communicator::getGrid2Rank( PartitionId pos[2], const PartitionId procgrid[2] ) const
{
    LAMA_ASSERT_ERROR( getSize() == procgrid[0] * procgrid[1], "illegal procgrid" )

    PartitionId rank = getRank();
    PartitionId size = procgrid[0];

    pos[1] = rank / size;
    rank = rank % size;
    pos[0] = rank;

    LAMA_LOG_INFO( logger,
                   *this << ": is (" << pos[0] << "," << pos[1] << ") of (" << procgrid[0] << "," << procgrid[1] << ")" )
}

void Communicator::getGrid3Rank( PartitionId pos[3], const PartitionId procgrid[3] ) const
{
    LAMA_ASSERT_ERROR( getSize() == procgrid[0] * procgrid[1] * procgrid[2], "illegal procgrid" )

    PartitionId rank = getRank();

    PartitionId size = procgrid[0] * procgrid[1];

    pos[2] = rank / size;
    rank = rank % size;

    size = procgrid[0];

    pos[1] = rank / size;
    rank = rank % size;
    pos[0] = rank;

    LAMA_LOG_INFO( logger,
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

        LAMA_LOG_INFO( logger,
                       "LAMA_NP=" << np4lama << " -> userProcArray " << userProcArray[0] << " x " << userProcArray[1] << " x " << userProcArray[2] )
    }
    else
    {
        LAMA_LOG_INFO( logger, "environment variable LAMA_NP no set" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename T>
auto_ptr<SyncToken> Communicator::defaultShiftAsync(
    T targetVals[],
    const T sourceVals[],
    const IndexType size,
    const int direction ) const
{
    // This default implementation uses synchronous shift and returns a NoSyncToken.

    IndexType recvSize = -1;
    recvSize = shiftImpl( targetVals, size, sourceVals, size, direction );

    LAMA_ASSERT_ERROR( recvSize == size, "asynchronous shift with different sizes on partitions" )

    return auto_ptr<SyncToken>( new NoSyncToken() );
}

std::auto_ptr<SyncToken> Communicator::shiftAsyncImpl(
    double newVals[],
    const double oldVals[],
    const IndexType size,
    const int direction ) const
{
    return defaultShiftAsync( newVals, oldVals, size, direction );
}

std::auto_ptr<SyncToken> Communicator::shiftAsyncImpl(
    float newVals[],
    const float oldVals[],
    const IndexType size,
    const int direction ) const
{
    return defaultShiftAsync( newVals, oldVals, size, direction );
}

std::auto_ptr<SyncToken> Communicator::shiftAsyncImpl(
    int newVals[],
    const int oldVals[],
    const IndexType size,
    const int direction ) const
{
    return defaultShiftAsync( newVals, oldVals, size, direction );
}

/* -------------------------------------------------------------------------- */

template<typename T>
IndexType Communicator::shift0(
    T targetVals[],
    const IndexType maxTargetSize,
    const T sourceVals[],
    const IndexType sourceSize ) const
{
    LAMA_ASSERT_ERROR( sourceSize <= maxTargetSize, "insufficient size for target array" )

    for ( IndexType i = 0; i < sourceSize; i++ )
    {
        targetVals[i] = sourceVals[i];
    }

    return sourceSize; // effective number of copied values.
}

/* -------------------------------------------------------------------------- */

template<typename T>
void Communicator::shiftArray( LAMAArray<T>& recvArray, const LAMAArray<T>& sendArray, const int direction ) const
{
    LAMA_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    if ( direction % getSize() == 0 )
    {
        // self assignment

        recvArray = sendArray;
        return;
    }

    HostReadAccess<T> sendData( sendArray );
    IndexType numSendElems = sendData.size();

    // make recv array large enough to fit for send data

    HostWriteOnlyAccess<T> recvData( recvArray, numSendElems );

    // but we are able to receive even more data if array is large enough

    IndexType maxNumRecvElems = recvData.capacity();

    // For shifting of data we use the pure virtual methods implemened by each communicator

    IndexType numRecvElems = shiftImpl( recvData.get(), maxNumRecvElems, sendData.get(), numSendElems, direction );

    LAMA_LOG_INFO( logger,
                   "shift, direction = " << direction << ", sent " << numSendElems << ", recvd " << numRecvElems << "( max was " << maxNumRecvElems << ")" )

    recvData.resize( numRecvElems ); // take over the size
}

/* -------------------------------------------------------------------------- */

template<typename T>
auto_ptr<SyncToken> Communicator::shiftAsync(
    LAMAArray<T>& recvArray,
    const LAMAArray<T>& sendArray,
    const int direction ) const
{
    LAMA_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    recvArray.clear(); // do not keep any old data, keep capacities

    auto_ptr<HostWriteAccess<T> > recvData( new HostWriteAccess<T>( recvArray ) );
    auto_ptr<HostReadAccess<T> > sendData( new HostReadAccess<T>( sendArray ) );

    IndexType numElems = sendData->size();

    recvData->resize( numElems ); // size should fit at least to keep own data

    // For shifting of data we use the pure virtual methods implemened by each communicator
    // Note: get is the method of the accesses and not of the auto_ptr

    auto_ptr<SyncToken> syncToken = shiftAsyncImpl( recvData->get(), sendData->get(), numElems, direction );

    LAMA_ASSERT_DEBUG( syncToken.get(), "NULL pointer for sync token" )

    // accesses are pushed in the sync token so they are freed after synchronization

    syncToken->pushAccess( auto_ptr<BaseAccess>( sendData.release() ) );
    syncToken->pushAccess( auto_ptr<BaseAccess>( recvData.release() ) );

    return syncToken;
}

/* -------------------------------------------------------------------------- */

template<typename T>
void Communicator::updateHalo( LAMAArray<T> &haloValues, const LAMAArray<T>& localValues, const Halo& halo ) const
{
    LAMA_REGION( "Communicator.updateHalo" )

    LAMA_LOG_INFO( logger, *this << ": update halo" )

    const CommunicationPlan& requiredPlan = halo.getRequiredPlan();

    LAMA_ASSERT_ERROR( requiredPlan.allocated(), "Required plan in Halo not allocated" )
    LAMA_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )

    const CommunicationPlan& providesPlan = halo.getProvidesPlan();

    LAMA_ASSERT_ERROR( providesPlan.allocated(), "Provides plan in Halo not allocated" )
    LAMA_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send

    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();

    LAMAArray<T> sendValues( numSendValues ); //!< temporary array for send communication

    LAMAArrayUtils::gather( sendValues, localValues, halo.getProvidesIndexes() );

    exchangeByPlan( haloValues, requiredPlan, sendValues, providesPlan );
}

/* -------------------------------------------------------------------------- */

template<typename T>
auto_ptr<SyncToken> Communicator::updateHaloAsync(
    LAMAArray<T>& haloValues,
    const LAMAArray<T>& localValues,
    const Halo& halo ) const
{
    LAMA_REGION( "Communicator.updateHaloAsync" )

    LAMA_LOG_INFO( logger, *this << ": asynchronous update halo" )

    const CommunicationPlan& requiredPlan = halo.getRequiredPlan();

    LAMA_ASSERT_ERROR( requiredPlan.allocated(), "Required plan in Halo not allocated" )
    LAMA_ASSERT_ERROR( requiredPlan.size() < getSize(), "Required plan in Halo mismatches size of communicator" )

    const CommunicationPlan& providesPlan = halo.getProvidesPlan();

    LAMA_ASSERT_ERROR( providesPlan.allocated(), "Provides plan in Halo not allocated" )
    LAMA_ASSERT_ERROR( providesPlan.size() < getSize(), "Provides plan in Halo mismatches size of communicator" )

    // Before we exchange by plan, we have to pack local values to send

    // Note: A previous MPI implementation took advantage of packing the data after
    //       starting the receives. This is here no more possible. But we might now
    //       pack the data already on the GPU and can avoid gpu->host transfer of all localValues

    IndexType numSendValues = providesPlan.totalQuantity();

    auto_ptr<LAMAArray<T> > sendValues( new LAMAArray<T>( numSendValues ) );

    // put together the (send) values to provide for other partitions

    {
        HostWriteAccess<T> sendData( *sendValues );
        HostReadAccess<T> localData( localValues );
        HostReadAccess<IndexType> sendIndexes( halo.getProvidesIndexes() );

        for ( IndexType i = 0; i < numSendValues; i++ )
        {
            sendData[i] = localData[sendIndexes[i]];
        }
    }

    auto_ptr<SyncToken> token = exchangeByPlanAsync( haloValues, requiredPlan, *sendValues, providesPlan );

    // Also push the sendValues array to the token so it will be freed after synchronization

    // Note: it is guaranteed that access to sendValues is freed before sendValues

    token->pushArray( auto_ptr<_LAMAArray>( sendValues.release() ) );

    return token;
}

/* -------------------------------------------------------------------------- */

void Communicator::computeOwners(
    const vector<IndexType>& requiredIndexes,
    const Distribution& distribution,
    vector<PartitionId>& owners ) const
{
    PartitionId rank = getRank();
    PartitionId size = getSize();

    const size_t requiredIndexesSize = requiredIndexes.size();

    LAMA_LOG_DEBUG( logger, "need owners for " << requiredIndexesSize << " global indexes" )

    if ( distribution.getCommunicator() != *this )
    {
        LAMA_THROWEXCEPTION( "The distribution has a different Communicator." )
    }

    int n = 0;
    owners.resize( requiredIndexesSize );

    //Check for own ownership. Mark needed Owners. Only exchange requests for unknown indexes.
    for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexesSize ); ++i )
    {
        if ( distribution.isLocal( requiredIndexes[i] ) )
        {
            owners[i] = rank;
        }
        else
        {
            n++;
            owners[i] = -1;
        }
    }

    LAMA_LOG_DEBUG( logger, requiredIndexesSize - n << " Indexes are local. Only need to send " << n << " values." )
    IndexType receiveSize = max( n ); // --> pure method call
    LAMA_LOG_DEBUG( logger, "max size of receive buffer is " << receiveSize )

    // Allocate the maxiamal needed size for the communication buffers

    LAMAArray<IndexType> indexesSendArray( receiveSize );
    LAMAArray<IndexType> indexesReceiveArray( receiveSize );
    LAMAArray<IndexType> ownersSendArray( receiveSize );
    LAMAArray<IndexType> ownersReceiveArray( receiveSize );
    {
        HostWriteAccess<IndexType> indexesSend( indexesSendArray );
        HostWriteAccess<IndexType> indexesReceive( indexesReceiveArray );
        HostWriteAccess<IndexType> ownersSend( ownersSendArray );
        HostWriteAccess<IndexType> ownersReceive( ownersReceiveArray );
        n = 0; // reset, counted again

        for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexesSize ); ++i )
        {
            if ( owners[i] == -1 )
            {
                indexesSend[n++] = requiredIndexes[i];
            }
        }

        // Important: set owners for full buffer of ownersSend

        for ( IndexType i = 0; i < receiveSize; ++i )
        {
            ownersSend[i] = -1;
        }
    }

    int ownersSize = -1;
    int currentSize = n;

    const int direction = 1; // send to right, recv from left
    for ( int iProc = 0; iProc < size - 1; ++iProc )
    {
        HostWriteAccess<IndexType> indexesSend( indexesSendArray );
        HostWriteAccess<IndexType> indexesReceive( indexesReceiveArray );
        HostWriteAccess<IndexType> ownersSend( ownersSendArray );
        HostWriteAccess<IndexType> ownersReceive( ownersReceiveArray );
        LAMA_LOG_DEBUG( logger,
                        *this << " shift: recv " << receiveSize << ", send " << currentSize << ", direction = " << direction )

        // --->   Pure method call

        currentSize = shiftImpl( indexesReceive.get(), receiveSize, indexesSend.get(), currentSize, direction );

        LAMA_ASSERT_ERROR( ownersSize == -1 || currentSize == ownersSize, "Communication corrupted." )

        LAMA_LOG_DEBUG( logger, "owners size = " << ownersSize << ", current size = " << currentSize )
        IndexType* indexes = indexesReceive.get();
        int* currentOwners = ownersSend.get();
        LAMA_LOG_DEBUG( logger, "check buffer with " << currentSize << " global indexes whether I am owner" )

        for ( int i = 0; i < currentSize; ++i )
        {
            //TODO there should be a blockwise implementation of isLocal
            LAMA_LOG_TRACE( logger,
                            "check global index " << indexes[i] << " with current owner " << currentOwners[i] << ", is local = " << distribution.isLocal( indexes[i] ) )

            if ( currentOwners[i] == -1 && distribution.isLocal( indexes[i] ) )
            {
                LAMA_LOG_TRACE( logger, *this << ": me is owner of global index " << indexes[i] )
                currentOwners[i] = rank;
            }
        }

        //Release so we can swap the Arrays
        indexesReceive.release();
        indexesSend.release();
        indexesReceiveArray.swap( indexesSendArray );

        LAMA_LOG_DEBUG( logger, *this << ": send array with " << currentSize << " owners to right" )

        for ( int i = 0; i < currentSize; i++ )
        {
            LAMA_LOG_TRACE( logger, *this << " send currentOwner[" << i << "] = " << ownersSend[i] )
        }

        // --->   Pure method call
        ownersSize = shiftImpl( ownersReceive.get(), receiveSize, ownersSend.get(), currentSize, direction );

        LAMA_LOG_DEBUG( logger, *this << ": recvd array with " << ownersSize << " owners from left" )
        for ( int i = 0; i < ownersSize; i++ )
        {
            LAMA_LOG_TRACE( logger, *this << ": recv currentOwner[" << i << "] = " << ownersReceive[i] )
        }

        //Release so we can swap the Arrays
        ownersSend.release();
        ownersReceive.release();
        ownersReceiveArray.swap( ownersSendArray );
    }

    HostWriteAccess<IndexType> ownersSend( ownersSendArray );
    for ( int i = 0; i < n; ++i )
    {
        LAMA_LOG_TRACE( logger,
                        *this << ": final " << i << " of " << n << ": " << requiredIndexes[i] << ", owner = " << ownersSend[i] )
    }

    // The Owner Indexes are always passed in the same order, so we can insert them easily.
    int nn = 0;

    for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexesSize ); ++i )
    {
        if ( owners[i] == -1 )
        {
            owners[i] = ownersSend[nn++];

            //TODO is this usefull for the speed ?
            if ( nn == n )
            {
                break;
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

// Instantiation of template methods for the supported types
template LAMA_DLL_IMPORTEXPORT
IndexType Communicator::shift0(
    float targetVals[],
    const IndexType maxTargetSize,
    const float sourceVals[],
    const IndexType sourceSize ) const;

template LAMA_DLL_IMPORTEXPORT
IndexType Communicator::shift0(
    double targetVals[],
    const IndexType maxTargetSize,
    const double sourceVals[],
    const IndexType sourceSize ) const;

template LAMA_DLL_IMPORTEXPORT
IndexType Communicator::shift0(
    int targetVals[],
    const IndexType maxTargetSize,
    const int sourceVals[],
    const IndexType sourceSize ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::shiftArray( LAMAArray<float>& recvArray, const LAMAArray<float>& sendArray, const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::shiftArray( LAMAArray<double>& recvArray, const LAMAArray<double>& sendArray, const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::shiftArray( LAMAArray<int>& recvArray, const LAMAArray<int>& sendArray, const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::defaultShiftAsync(
    double targetVals[],
    const double sourceVals[],
    const IndexType size,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::defaultShiftAsync(
    float targetVals[],
    const float sourceVals[],
    const IndexType size,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::defaultShiftAsync(
    int targetVals[],
    const int sourceVals[],
    const IndexType size,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::shiftAsync(
    LAMAArray<float>& recvArray,
    const LAMAArray<float>& sendArray,
    const int direction ) const;
template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::shiftAsync(
    LAMAArray<double>& recvArray,
    const LAMAArray<double>& sendArray,
    const int direction ) const;
template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::shiftAsync(
    LAMAArray<int>& recvArray,
    const LAMAArray<int>& sendArray,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::updateHalo(
    LAMAArray<float>& haloValues,
    const LAMAArray<float>& localValues,
    const Halo& halo ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::updateHalo(
    LAMAArray<double>& haloValues,
    const LAMAArray<double>& localValues,
    const Halo& halo ) const;

template LAMA_DLL_IMPORTEXPORT
void Communicator::updateHalo( LAMAArray<int>& haloValues, const LAMAArray<int>& localValues, const Halo& halo ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::updateHaloAsync(
    LAMAArray<float>& haloValues,
    const LAMAArray<float>& localValues,
    const Halo& halo ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::updateHaloAsync(
    LAMAArray<double>& haloValues,
    const LAMAArray<double>& localValues,
    const Halo& halo ) const;

template LAMA_DLL_IMPORTEXPORT
auto_ptr<SyncToken> Communicator::updateHaloAsync(
    LAMAArray<int>& haloValues,
    const LAMAArray<int>& localValues,
    const Halo& halo ) const;

}
