/**
 * @file Distribution.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of methods for base class Distribution.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/dmemo/Distribution.hpp>

// local library
#include <scai/dmemo/Distributed.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/typeloop.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

using namespace hmemo;
using common::TypeTraits;

namespace dmemo
{

/* ------  Static class variables --------------------------------------- */

SCAI_LOG_DEF_LOGGER( Distribution::logger, "Distribution" )

/* ------  Constructor  ------------------------------------------------- */

Distribution::Distribution( const IndexType globalSize )
    : mGlobalSize( globalSize ), mCommunicator( Communicator::getCommunicatorPtr( Communicator::NO ) )
{
    SCAI_LOG_INFO( logger, "Distribution(" << mGlobalSize << ") onto NoCommunicator" )
}

/* ------  Constructor  ------------------------------------------------- */

Distribution::Distribution( const IndexType globalSize, const CommunicatorPtr communicator )
    : mGlobalSize( globalSize ), mCommunicator( communicator )
{
    if ( !mCommunicator )
    {
        COMMON_THROWEXCEPTION( "Distribution without a Communicator is not allowed" )
    }

    SCAI_LOG_INFO( logger, "Distribution( size = " << globalSize << ", comm = " << *mCommunicator << " )" )
}

/* -------  Destructor  ------------------------------------------------- */

Distribution::~Distribution()
{
    SCAI_LOG_DEBUG( logger, "~Distribution, global size " << mGlobalSize << ", communicator = " << mCommunicator )
}

/* ---------------------------------------------------------------------- */

bool Distribution::operator==( const Distribution& other ) const
{
    // two distributions are the same if they are both replicated
    SCAI_LOG_TRACE( logger, "check " << *this << " == " << other )
    bool isSame = false;

    if ( this == &other )
    {
        isSame = true;
        SCAI_LOG_DEBUG( logger, *this << " == " << other << ": pointer equal" )
    }
    else if ( isReplicated() && other.isReplicated() )
    {
        isSame = getGlobalSize() == other.getGlobalSize();
        SCAI_LOG_DEBUG( logger, *this << " == " << other << ": both are replicated, same size" )
    }
    else
    {
        isSame = isEqual( other );
        SCAI_LOG_DEBUG( logger, *this << " == " << other << ": " << isSame )
    }

    return isSame;
}

/* ---------------------------------------------------------------------- */

bool Distribution::operator!=( const Distribution& other ) const
{
    return !( *this == other );
}

/* ---------------------------------------------------------------------- */

const Communicator& Distribution::getCommunicator() const
{
    return *mCommunicator;
}

/* ---------------------------------------------------------------------- */

CommunicatorPtr Distribution::getCommunicatorPtr() const
{
    return mCommunicator;
}

/* ---------------------------------------------------------------------- */

PartitionId Distribution::getNumPartitions() const
{
    // mCommunicator is never NULL, but just in case
    SCAI_ASSERT( mCommunicator, "Distribution without a Communicator is not allowed" )
    return mCommunicator->getSize();
}

/* ---------------------------------------------------------------------- */

void Distribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "Distribution";
}

/* ---------------------------------------------------------------------- */

void Distribution::computeOwners(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    SCAI_LOG_INFO( logger, "compute owners via communicator (default)" )
    // use communicator to compute ownership on each processor
    IndexType n = requiredIndexes.size();
    owners.resize( n );
    mCommunicator->computeOwners( &owners[0], *this, &requiredIndexes[0], n );
}

/* ---------------------------------------------------------------------- */

template<typename T1, typename T2>
void Distribution::replicate( T1* allValues, const T2* localValues ) const
{
    SCAI_REGION( "Distribution.replicate" )

    const Communicator& comm = getCommunicator();
    IndexType currentSize = getLocalSize();
    // Implemenation via cyclic shifting of the vector data and distribution
    IndexType maxLocalSize = comm.max( currentSize );

    SCAI_LOG_INFO( logger,
                   comm << ": replicate localValues<" << TypeTraits<T2>::id() << ">[ " << currentSize
                   << ", max = " << maxLocalSize << " ] " << " to allValues<" << TypeTraits<T1>::id() << ">[ " << getGlobalSize() << " ]" )

    // Only allocate the needed size of the Arrays

    HArray<T1> valuesSend;
    HArray<T1> valuesReceive;
    HArray<IndexType> indexesSend;
    HArray<IndexType> indexesReceive;

    IndexType countValues = 0; // count values set in the global vector

    // set my owned indexes and my values

    ContextPtr commContext = comm.getCommunicationContext( valuesSend );

    // capacity of send arrays should also be sufficient for receiving data

    indexesSend.reserve( commContext, maxLocalSize );
    valuesSend.reserve( commContext, maxLocalSize );

    SCAI_LOG_INFO( logger, "replicate on this communication context: " << *commContext )

    {
        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentSize );
        WriteOnlyAccess<T1> wValuesSend( valuesSend, commContext, currentSize );

        // current workaround as commContext works like HostContext

        IndexType* pIndexesSend = wIndexesSend.get();
        T1* pValuesSend = wValuesSend.get();

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            IndexType globalIndex = local2global( i );
            SCAI_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )
            pIndexesSend[i] = globalIndex;
            pValuesSend[i] = static_cast<T1>( localValues[i] ); // type conversion here
            allValues[globalIndex] = static_cast<T1>( localValues[i] ); // type conversion here
        }
    }

    // capacity of receive arrays should be sufficient for receiving data to avoid reallocations

    indexesReceive.reserve( commContext, maxLocalSize );
    valuesReceive.reserve( commContext, maxLocalSize );

    countValues += currentSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1 );
        comm.shiftArray( valuesReceive, valuesSend, 1 );

        SCAI_ASSERT_EQ_DEBUG( valuesReceive.size(), indexesReceive.size(), "size mismatch" )

        currentSize = valuesReceive.size();

        // sort in the received values

        {
            ReadAccess<IndexType> readIndexesReceive( indexesReceive, commContext );
            ReadAccess<T1> readValuesReceive( valuesReceive, commContext );

            // current workaround as commContext works like HostContext

            const IndexType* rIndexesReceive = readIndexesReceive.get();
            const T1* rValuesReceive = readValuesReceive.get();

            for ( IndexType i = 0; i < currentSize; i++ )
            {
                const IndexType globalIndex = rIndexesReceive[i];
                SCAI_ASSERT_DEBUG( globalIndex < getGlobalSize(),
                                   *this << ": global index " << globalIndex << " illegal" )
                allValues[globalIndex] = rValuesReceive[i]; // implicit type conversion done here
            }
        }

        countValues += currentSize;
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    SCAI_LOG_INFO( logger, "replicated by " << ( np - 1 ) << " array shifts" )

    // # globalSize values must have been counted
    SCAI_ASSERT_EQ_DEBUG( countValues, getGlobalSize(), "" )
}

/* ---------------------------------------------------------------------- */

template<typename T1, typename T2>
void Distribution::replicateN( T1* allValues, const T2* localValues, const IndexType n ) const
{
    SCAI_REGION( "Distribution.replicateN" )

    const Communicator& comm = getCommunicator();

    // Implemenation via cyclic shifting of the vector data and distribution

    // maximal number of elements needed to allocate sufficient receive buffer but also to avoid reallocations

    IndexType currentSize = getLocalSize();
    IndexType maxLocalSize = comm.max( currentSize );

    SCAI_LOG_INFO( logger,
                   comm << ": replicateN, n = " << n << ", localValues<" << common::getScalarType<T2>() << ">[ " << currentSize << ", max = " << maxLocalSize << " ] " << " to allValues<" << common::getScalarType<T1>() << ">[ " << getGlobalSize() << " ]" )

    // Only allocate the needed size of the Arrays

    HArray<T1> valuesSend;
    HArray<T1> valuesReceive;
    HArray<IndexType> indexesSend;
    HArray<IndexType> indexesReceive;

    // set my owned indexes and my values

    ContextPtr commContext = comm.getCommunicationContext( valuesSend );

    // capacity of send arrays should also be sufficient for receiving data

    indexesSend.reserve( commContext, maxLocalSize );
    valuesSend.reserve( commContext, maxLocalSize * n );

    SCAI_LOG_INFO( logger, "replicate on this communication context: " << *commContext )

    {
        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentSize );
        WriteOnlyAccess<T1> wValuesSend( valuesSend, commContext, currentSize * n );

        // current workaround as commContext works like HostContext

        IndexType* pIndexesSend = wIndexesSend.get();
        T1* pValuesSend = wValuesSend.get();

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            IndexType globalIndex = local2global( i );
            SCAI_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )
            pIndexesSend[i] = globalIndex;

            for ( IndexType j = 0; j < n; j++ )
            {
                pValuesSend[i * n + j] = static_cast<T1>( localValues[i * n + j] ); // type conversion here
                allValues[globalIndex * n + j] = static_cast<T1>( localValues[i * n + j] ); // type conversion here
            }

            pValuesSend[i] = static_cast<T1>( localValues[i] ); // type conversion here
            allValues[globalIndex] = static_cast<T1>( localValues[i] ); // type conversion here
        }
    }

    IndexType countLines = currentSize; // count values set in the global vector, currentSize has been done

    // capacity of receive arrays should be sufficient for receiving data to avoid reallocations

    indexesReceive.reserve( commContext, maxLocalSize );
    valuesReceive.reserve( commContext, maxLocalSize * n );

    // now nproc - 1 steps for cyclic shifting

    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1 );
        currentSize = indexesReceive.size(); // next numbe of lines to deal with

        comm.shiftArray( valuesReceive, valuesSend, 1 );

        SCAI_ASSERT_EQ_DEBUG( currentSize * n, valuesReceive.size(), "size mismatch for replicate n = " << n )

        // sort in the received values

        {
            ReadAccess<IndexType> readIndexesReceive( indexesReceive, commContext );
            ReadAccess<T1> readValuesReceive( valuesReceive, commContext );

            // current workaround as commContext works like HostContext

            const IndexType* rIndexesReceive = readIndexesReceive.get();
            const T1* rValuesReceive = readValuesReceive.get();

            for ( IndexType i = 0; i < currentSize; i++ )
            {
                const IndexType globalIndex = rIndexesReceive[i];
                SCAI_ASSERT_DEBUG( globalIndex < getGlobalSize(),
                                   *this << ": global index " << globalIndex << " illegal" )

                for ( IndexType j = 0; j < n; j++ )
                {
                    allValues[globalIndex * n + j] = rValuesReceive[i * n + j];
                }
            }
        }

        countLines += currentSize;
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    SCAI_ASSERT_EQ_DEBUG( countLines, getGlobalSize(), "not all lines seen" )
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
static IndexType fillGlobal(
    ValueType* allValues,
    const IndexType* allOffsets,
    const IndexType* indexes,
    const IndexType numIndexes,
    const ValueType* values )
{
    IndexType counter = 0; // traverses values, counts the number of filled values

    for ( IndexType i = 0; i < numIndexes; i++ )
    {
        const IndexType allIndex = indexes[i];
        const IndexType offset = allOffsets[allIndex];
        const IndexType size = allOffsets[allIndex + 1] - offset;

        /* No logger available here:

         std::cout << "fill index " << i << " of " << numIndexes << ": global = "
         << allIndex << ", size = " << size << ", offset = " << offset  << std::endl;
         */

        for ( IndexType j = 0; j < size; j++ )
        {
            allValues[offset + j] = values[counter + j];
        }

        counter += size;
    }

    return counter; // number of filled elements
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
void Distribution::replicateRagged(
    ValueType allValues[],
    const ValueType localValues[],
    const IndexType allOffsets[] ) const
{
    IndexType currentElemSize = getLocalSize();
    const Communicator& comm = getCommunicator();

    // we need the offsets for allValues to store directly the values
    IndexType maxLocalSize = comm.max( currentElemSize );

    HArray<IndexType> indexesSend;
    HArray<IndexType> indexesReceive;

    ContextPtr commContext = comm.getCommunicationContext( indexesSend );

    indexesReceive.reserve( commContext, maxLocalSize );
    indexesSend.reserve( commContext, maxLocalSize );

    IndexType currentDataSize = 0;

    {
        // due to currentElemSize <= maxLocalSize

        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentElemSize );

        IndexType* pIndexesSend = wIndexesSend.get(); // is a HostContext

        // pack my indexes, work as global indexes for allValues

        for ( IndexType i = 0; i < currentElemSize; i++ )
        {
            pIndexesSend[i] = local2global( i );
        }

        // fill my local values in all values
        currentDataSize = fillGlobal( allValues, allOffsets, pIndexesSend, currentElemSize, localValues );
    }

    SCAI_LOG_DEBUG( logger,
                    comm << ": filled my local data: " << currentElemSize << " buckets with " << currentDataSize << " values" )
    // get maximal size of data values and allocate send buffer
    IndexType maxLocalDataSize = comm.max( currentDataSize );
    // Implemenation via cyclic shifting of the vector data and distribution
    SCAI_LOG_DEBUG( logger, "maximal data size for exchange = " << maxLocalDataSize )

    SCAI_LOG_INFO( logger,
                   comm << ": replicateRagged<" << common::getScalarType<ValueType>() << ">, localValues[ " << currentDataSize << ", max = " << maxLocalDataSize << " ] " << " to allValues [ allOffsets[ " << getGlobalSize() << " ] ]" )

    HArray<ValueType> valuesSend;
    HArray<ValueType> valuesReceive;

    valuesSend.reserve( commContext, maxLocalDataSize );
    valuesReceive.reserve( commContext, maxLocalDataSize );

    {
        WriteOnlyAccess<ValueType> wValuesSend( valuesSend, commContext, currentDataSize );

        ValueType* pValuesSend = wValuesSend.get();

        // fill my local values in send buffer

        for ( IndexType i = 0; i < currentDataSize; i++ )
        {
            pValuesSend[i] = localValues[i];
        }
    }

    // Define counters for exchanged values
    IndexType countElemValues = currentElemSize;
    IndexType countDataValues = currentDataSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1 );
        comm.shiftArray( valuesReceive, valuesSend, 1 );

        IndexType newSize1 = indexesReceive.size();
        IndexType newSize2 = valuesReceive.size();

        // copy the received values in allValues at the right place

        {
            ReadAccess<IndexType> rIndexesReceive( indexesReceive, commContext );
            ReadAccess<ValueType> rValuesReceive( valuesReceive, commContext );

            IndexType size = -1;
            size = fillGlobal( allValues, allOffsets, rIndexesReceive.get(), newSize1, rValuesReceive.get() );
            SCAI_LOG_DEBUG( logger,
                            comm << ": filled received data: " << newSize1 << " buckets with " << size << " values" )
            SCAI_ASSERT_EQ_DEBUG( size, newSize2, "size mismatch" )
        }

        currentElemSize = newSize1;
        currentDataSize = newSize2;
        countElemValues += currentElemSize;
        countDataValues += currentDataSize;
        // swap the send and receive buffer
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    // verify that all values are available

    SCAI_ASSERT_EQ_DEBUG( countElemValues, getGlobalSize(), "not all elems seen" )
    SCAI_ASSERT_EQ_DEBUG( countDataValues, allOffsets[getGlobalSize()], "not all data seen" )
}

/* ---------------------------------------------------------------------- */

Distribution* Distribution::getDistribution(
    const std::string& kind,
    CommunicatorPtr comm,
    const IndexType globalSize,
    const float weight )
{
    return Distribution::create( kind, DistributionArguments( comm, globalSize, NULL, weight ) );
}

Distribution* Distribution::getDistribution(
    const std::string& kind,
    CommunicatorPtr comm,
    const Distributed& matrix,
    const float weight )
{
    IndexType globalSize = matrix.getDistribution().getGlobalSize();

    return Distribution::create( kind, DistributionArguments( comm, globalSize, &matrix, weight ) );
}

/* ---------------------------------------------------------------------- */

#define DMEMO_DISTRIBUTE2_INST( ValueType, OtherValueType )                                                                         \
    template void Distribution::replicateN<ValueType, OtherValueType>( ValueType*, const OtherValueType*, const IndexType ) const;  \
    template void Distribution::replicate<ValueType, OtherValueType>( ValueType*, const OtherValueType* ) const;

#define DMEMO_DISTRIBUTE_INST( ValueType )  \
    template void Distribution::replicateRagged<ValueType>( ValueType*, const ValueType*, const IndexType* ) const;                 \
    SCAI_COMMON_TYPELOOP_LVL2( ValueType, DMEMO_DISTRIBUTE2_INST, SCAI_ARITHMETIC_ARRAY_HOST )

SCAI_COMMON_TYPELOOP( DMEMO_DISTRIBUTE_INST, SCAI_ARITHMETIC_ARRAY_HOST )
     
// template instantiation for the supported data types

#undef DMEMO_DISTRIBUTE2_INST
#undef DMEMO_DISTRIBUTE_INST
/* ---------------------------------------------------------------------- */

} /* end namespace dmemo */

} /* end namespace scai */
