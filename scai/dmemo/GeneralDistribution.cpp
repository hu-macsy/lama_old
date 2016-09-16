/**
 * @file GeneralDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief GeneralDistribution.cpp
 * @author brandes
 * @date 25.02.2011
 */

// hpp
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>

// internal scai libraries
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/common/macros/assert.hpp>

#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

// std
#include <algorithm>
#include <functional>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GeneralDistribution::logger, "Distribution.General" )

const char GeneralDistribution::theCreateValue[] = "GENERAL";

GeneralDistribution::GeneralDistribution(
    const IndexType globalSize,
    const HArray<IndexType>& myIndexes,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator )
{
    ReadAccess<IndexType> rIndexes( myIndexes );

    IndexType nLocal = myIndexes.size();

    WriteOnlyAccess<IndexType> wLocal2Global( mLocal2Global, nLocal );

    for ( IndexType localIndex = 0; localIndex < nLocal; ++localIndex )
    {
        IndexType globalIndex = rIndexes[ localIndex ];

        SCAI_ASSERT_LT_ERROR( globalIndex, mGlobalSize, "global index out of range" )

        wLocal2Global[ localIndex ]  = globalIndex;
        mGlobal2Local[ globalIndex ] = localIndex;
    }

    // Note: the constructor is completely local, but make some consistency check now 

    SCAI_ASSERT_EQ_ERROR( mGlobalSize, communicator->sum( nLocal ), "illegal general distribution" )
}

GeneralDistribution::GeneralDistribution(
    const HArray<PartitionId>& owners,
    const CommunicatorPtr communicator ) : 

    Distribution( 0, communicator )

{
    using namespace hmemo;

    ContextPtr ctx = Context::getHostPtr();

    PartitionId rank = mCommunicator->getRank();
    PartitionId size = mCommunicator->getSize();

    if ( rank == MASTER )
    {
        mGlobalSize = owners.size();
    }

    mCommunicator->bcast( &mGlobalSize, 1, MASTER );

    SCAI_LOG_DEBUG( logger, *mCommunicator << ": global size = " << mGlobalSize )

    utilskernel::LArray<IndexType> localSizes;
    utilskernel::LArray<IndexType> localOffsets;

    // count in localSizes for each partition the owners
    // owners = [ 0, 1, 2, 1, 2, 1, 0 ] -> sizes = [2, 3, 2 ]
    // sizes.sum() == owners.size()

    if ( rank == MASTER )
    {
        utilskernel::HArrayUtils::bucketCount( localSizes, owners, size );
        IndexType lsum = localSizes.sum();
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": sum( localSizes ) = " << lsum << ", must be " << mGlobalSize );
    }
    else
    { 
        // all other procs intialize localSizes with a dummy value to avoid read access of uninitialized array
        utilskernel::HArrayUtils::setOrder( localSizes, 1 );
    }

    IndexType localSize;  

    // scatter partition sizes

    {
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": before scatter, localSizes = " << localSizes )
        ReadAccess<IndexType> rSizes( localSizes );
        mCommunicator->scatter( &localSize, 1, MASTER, rSizes.get() );
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": after scatter, localSize = " << localSize )
    }

    SCAI_LOG_DEBUG( logger, *mCommunicator << ": owns " << localSize << " of " << mGlobalSize << " elements" )

    if ( rank == MASTER )
    {
        localOffsets.reserve( ctx, size + 1 );
        localOffsets = localSizes;
        utilskernel::HArrayUtils::scan( localOffsets );
        SCAI_LOG_DEBUG( logger, "scan done, sum = " << localOffsets[ size ] )
    }

    // Now resort 0, ..., n - 1 according to the owners

    HArray<IndexType> sortedIndexes;

    if ( rank == MASTER )
    {
        SCAI_LOG_DEBUG( logger, "reorder for indexes" )

        ContextPtr loc = Context::getHostPtr();

        static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::sortInBuckets<PartitionId> > sortInBuckets;

        sortInBuckets.getSupportedContext( loc );

        WriteOnlyAccess<IndexType> wIndexes( sortedIndexes, loc, mGlobalSize );
        WriteAccess<IndexType> wOffsets( localOffsets, loc );
        ReadAccess<PartitionId> rOwners( owners, loc );

        sortInBuckets[loc]( wIndexes, wOffsets, size, rOwners, mGlobalSize );
    }
    else
    { 
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": initialize sortedIndexes " )

        // all other procs intialize sortedIndexes with a dummy value to avoid read access of uninitialized array
        utilskernel::HArrayUtils::setOrder( sortedIndexes, 1 );
    }

    WriteOnlyAccess<IndexType> wLocal2Global( mLocal2Global, localSize );

    {
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": before scatterV, sortedIndexes = " << sortedIndexes  )
        ReadAccess<IndexType> rIndexes( sortedIndexes );
        ReadAccess<IndexType> rSizes( localSizes );
        mCommunicator->scatterV( wLocal2Global.get(), localSize, MASTER, rIndexes.get(), rSizes.get() );
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": after scatterV, sortedIndexes = " << sortedIndexes )
    }

    // Compute Global2Local

    for ( IndexType i = 0; i < localSize; ++i )
    {
        mGlobal2Local[ wLocal2Global[i] ] = i;
    }
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution( const Distribution& other ) : 

    Distribution( other.getGlobalSize(), other.getCommunicatorPtr() )

{
    WriteOnlyAccess<IndexType> wLocal2Global( mLocal2Global, other.getLocalSize() );

    for ( IndexType i = 0; i < getGlobalSize(); ++i )
    {
        if ( other.isLocal( i ) )
        {
            IndexType localIndex = other.global2local( i );
            mGlobal2Local[i] = localIndex;
            wLocal2Global[localIndex] = i;
        }
    }
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution( const GeneralDistribution& other ) :

    Distribution( other.getGlobalSize(), other.getCommunicatorPtr() )

{
    mLocal2Global = other.mLocal2Global;
    mGlobal2Local = other.mGlobal2Local;
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution( const IndexType globalSize, const CommunicatorPtr communicator ) : 

    Distribution( globalSize, communicator )
{
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::~GeneralDistribution()
{
    SCAI_LOG_INFO( logger, "~GeneralDistribution" )
}

bool GeneralDistribution::isLocal( const IndexType index ) const
{
    return mGlobal2Local.find( index ) != mGlobal2Local.end();
}

IndexType GeneralDistribution::getLocalSize() const
{
    return static_cast<IndexType>( mLocal2Global.size() );
}

IndexType GeneralDistribution::local2global( const IndexType localIndex ) const
{
    return mLocal2Global[localIndex];
}

IndexType GeneralDistribution::global2local( const IndexType globalIndex ) const
{
    const Global2LocalMapType::const_iterator elem = mGlobal2Local.find( globalIndex );

    if ( elem == mGlobal2Local.end() )
    {
        return nIndex;
    }

    return elem->second;
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getBlockDistributionSize() const
{
    // Note: we assume that the local indexes are descending and do not contain doubles
 
    IndexType localSize = mLocal2Global.size();

    bool isBlocked = true;

    ReadAccess<IndexType> rIndexes( mLocal2Global );

    for ( IndexType i = 0; i < localSize; ++i )
    {
        if ( rIndexes[i] - rIndexes[0] != i )
        {
            isBlocked = false;
        }
    }

    CommunicatorPtr comm = getCommunicatorPtr();

    isBlocked = comm->all( isBlocked );

    if ( !isBlocked )
    {
        return nIndex;
    }

    // Each processor has a contiguous part, but verify that it is in the same order

    GenBlockDistribution genBlock( mGlobalSize, localSize, comm );

    IndexType lb;
    IndexType ub;

    genBlock.getLocalRange( lb, ub );

    isBlocked = true;

    if ( localSize > 0 )
    {
        isBlocked = ( rIndexes[0] == lb ) && ( rIndexes[localSize - 1] == ub );
    }

    isBlocked = comm->all( isBlocked );

    if ( !isBlocked )
    {
        return nIndex;
    }

    return localSize;
}

/* ---------------------------------------------------------------------- */

bool GeneralDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( other.getKind() != getKind() )
    {
        return false;
    }

    const GeneralDistribution& otherGen = reinterpret_cast<const GeneralDistribution&>( other );

    bool localSameSize = otherGen.getLocalSize() == getLocalSize();
    bool allSameSize = mCommunicator->all( localSameSize );

    SCAI_LOG_DEBUG( logger, *this << ": localSameSize = " << localSameSize << ", allSameSize = " << allSameSize )
 
    if ( !allSameSize ) 
    {
        return false;
    }

    // values will only be compared it sizes are same on all processors

    bool localSameVals = mLocal2Global.maxDiffNorm( otherGen.getMyIndexes() ) == 0;
    bool allSameVals   = mCommunicator->all( localSameVals );

    SCAI_LOG_DEBUG( logger, *this << ": localSameVals = " << localSameVals << ", allSameVals = " << allSameVals )

    return allSameVals;
}

void GeneralDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "GeneralDistribution( size = " << mLocal2Global.size() << " of " << mGlobalSize << ", comm = "
           << *mCommunicator << " )";
}

static void setOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes, const HArray<IndexType>& offsets )
{
    IndexType globalSize = indexes.size();
    IndexType nOwners    = offsets.size() - 1;

    WriteOnlyAccess<PartitionId> wOwners( owners, indexes.size() );
    ReadAccess<IndexType> rOffsets( offsets );
    ReadAccess<IndexType> rIndexes( indexes );

    // for testing: init

    for ( IndexType i = 0; i < globalSize; ++i )
    {
        wOwners[i] = nPartition;
    }

    for ( IndexType owner = 0; owner < nOwners; ++owner )
    {
        for ( IndexType j = rOffsets[owner]; j < rOffsets[owner + 1]; ++j )
        {
            wOwners[rIndexes[j]] = owner;
        }
    }
}

void GeneralDistribution::allOwners( HArray<PartitionId>& owners, const PartitionId root ) const
{
    ContextPtr ctx = Context::getHostPtr();

    PartitionId rank  = mCommunicator->getRank();
    PartitionId parts = mCommunicator->getSize();

    IndexType localSize = getLocalSize();

    // gather number of local rows

    HArray<IndexType> localSizes;

    {
        WriteOnlyAccess<IndexType> wLocalSizes( localSizes, ctx, parts );
        mCommunicator->gather( wLocalSizes.get(), 1, root, &localSize );
    }

    HArray<IndexType> offsets;
    offsets.reserve( ctx, parts + 1 );;

    if ( rank == root )
    {
        utilskernel::HArrayUtils::assign( offsets, localSizes, ctx );
        IndexType nTotal = utilskernel::HArrayUtils::scan( offsets );
        SCAI_ASSERT( nTotal == mGlobalSize, "sum of local rows is not global size" )
    }

    // gather global indices of local rows

    HArray<IndexType> allIndexes( rank == root ? mGlobalSize : 2 );
   
    {   
        WriteAccess<IndexType> wAllIndexes( allIndexes );
        ReadAccess<IndexType> rLocalSizes( localSizes );
        ReadAccess<IndexType> rIndexes( mLocal2Global );
        mCommunicator->gatherV( wAllIndexes.get(), localSize, root, rIndexes.get(), rLocalSizes.get() );
    }

    // now we ca scatter the owner ids in the owner array

    if ( rank == root )
    {
        setOwners( owners, allIndexes, offsets );
    }
}

void GeneralDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    utilskernel::HArrayUtils::assign( myGlobalIndexes, mLocal2Global );
}

} /* end namespace dmemo */

} /* end namespace scai */
