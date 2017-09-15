/**
 * @file Redistributor.cpp
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
 * @brief Implementation of methods for Redistributor class.
 * @author Thomas Brandes
 * @date 08.10.2011
 */

// hpp
#include <scai/dmemo/Redistributor.hpp>

// local library
#include <scai/dmemo/HaloBuilder.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/unique_ptr.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::unique_ptr;
using common::scoped_array;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( Redistributor::logger, "Redistributor" )

Redistributor::Redistributor( DistributionPtr targetDistribution, DistributionPtr sourceDistribution )

    : mSourceDistribution( sourceDistribution ), mTargetDistribution( targetDistribution )

{
    SCAI_ASSERT_ERROR( sourceDistribution, "NULL pointer for source distribution" )
    SCAI_ASSERT_ERROR( targetDistribution, "NULL pointer for target distribution" )
    // Dereference distribution pointers, avoids checks for validity
    const Distribution& sourceDist = *sourceDistribution;
    const Distribution& targetDist = *targetDistribution;
    SCAI_ASSERT_EQ_ERROR( sourceDist.getGlobalSize(), targetDist.getGlobalSize(), "serious size mismatch" )
    SCAI_ASSERT_EQ_ERROR( sourceDist.getCommunicator(), targetDist.getCommunicator(), "redistribute only with same communicator" )
    mSourceSize = sourceDist.getLocalSize();
    mTargetSize = targetDist.getLocalSize();
    SCAI_LOG_INFO( logger,
                   sourceDist.getCommunicator() << ": build redistributor " << targetDist << " <- " << sourceDist << ", have " << mSourceSize << " source values and " << mTargetSize << " target values" )
    // localSourceIndexes, localTargetIndexes are used for local permutation
    // we do not know the exact sizes now so we take maximal value
    WriteOnlyAccess<IndexType> localSourceIndexes( mLocalSourceIndexes, mSourceSize );
    WriteOnlyAccess<IndexType> localTargetIndexes( mLocalTargetIndexes, mTargetSize );
    std::vector<IndexType> requiredIndexes;
    mNumLocalValues = 0; // count number of local copies from source to target

    for ( IndexType i = 0; i < mTargetSize; i++ )
    {
        IndexType globalIndex = targetDist.local2global( i );

        if ( sourceDist.isLocal( globalIndex ) )
        {
            IndexType sourceLocalIndex = sourceDist.global2local( globalIndex );
            SCAI_LOG_TRACE( logger,
                            "target local index " << i << " is global " << globalIndex << ", is source local index " << sourceLocalIndex )
            //  so globalIndex is local in both distributions
            localTargetIndexes[mNumLocalValues] = i; // where to scatter in target
            localSourceIndexes[mNumLocalValues] = sourceLocalIndex;
            mNumLocalValues++;
        }
        else
        {
            SCAI_LOG_TRACE( logger, "target local index " << i << " is global " << globalIndex << ", is remote" )
            // needed from other processor
            requiredIndexes.push_back( globalIndex );
        }
    }

    // Adapt sizes of arrays with local indexes
    localSourceIndexes.resize( mNumLocalValues );
    localTargetIndexes.resize( mNumLocalValues );
    SCAI_LOG_DEBUG( logger,
                    sourceDist.getCommunicator() << ": target dist has local " << mTargetSize << " vals, " << mNumLocalValues << " are local, " << requiredIndexes.size() << " are remote." )
    // Halo is only for exchange of non-local values
    HArrayRef<IndexType> arrRequiredIndexes( requiredIndexes );
    HaloBuilder::build( sourceDist, arrRequiredIndexes, mHalo );
    // Set in the source index vector the values to provide for other processors
    const Halo& halo = mHalo;
    SCAI_LOG_INFO( logger,
                   sourceDist.getCommunicator() << ": halo source has " << halo.getProvidesPlan().totalQuantity() << " indexes, " << "halo target has " << halo.getRequiredPlan().totalQuantity() << " indexes" )
    const CommunicationPlan& providesPlan = halo.getProvidesPlan();
    ContextPtr contextPtr = Context::getHostPtr();
    WriteAccess<IndexType> haloSourceIndexes( mHaloSourceIndexes, contextPtr );
    WriteAccess<IndexType> haloTargetIndexes( mHaloTargetIndexes, contextPtr );
    ReadAccess<IndexType> haloProvidesIndexes( halo.getProvidesIndexes(), contextPtr );
    haloSourceIndexes.resize( providesPlan.totalQuantity() );
    IndexType offset = 0; // runs through halo source indexes

    for ( PartitionId i = 0; i < providesPlan.size(); i++ )
    {
        IndexType n = providesPlan[i].quantity;

        const IndexType* pindexes = haloProvidesIndexes.get() + providesPlan[i].offset;

        for ( IndexType j = 0; j < n; j++ )
        {
            SCAI_LOG_TRACE( logger, "halo source index[" << offset << "] = " << pindexes[j] )
            haloSourceIndexes[offset++] = pindexes[j];
        }
    }

    SCAI_LOG_INFO( logger, "have set " << offset << " halo source indexes" )
    // In contrary to Halo schedules we have here the situation that each non-local
    // index of source should be required by some other processor.
    SCAI_ASSERT_EQ_ERROR( offset, haloSourceIndexes.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mNumLocalValues + offset, mSourceSize, "serious mismatch" )
    // Now add the indexes where to scatter the halo into destination
    IndexType haloSize = halo.getHaloSize();
    SCAI_ASSERT_ERROR( mNumLocalValues + haloSize == mTargetSize, "size mismatch" )
    haloTargetIndexes.resize( haloSize );

    for ( IndexType i = 0; i < haloSize; i++ )
    {
        IndexType globalIndex = requiredIndexes[i];
        IndexType localIndex = targetDist.global2local( globalIndex );
        IndexType haloIndex = mHalo.global2halo( globalIndex );
        SCAI_LOG_TRACE( logger,
                        "saved mapping for target dist: local = " << localIndex << ", global = " << globalIndex << ", halo = " << haloIndex )
        haloTargetIndexes[haloIndex] = localIndex;
    }

    // all info needed for redistribution is now available
    SCAI_LOG_INFO( logger, "constructed " << *this )
}

Redistributor::Redistributor( const HArray<IndexType>& newOwners, DistributionPtr sourceDistribution )

    : mSourceDistribution( sourceDistribution )
{
    SCAI_ASSERT_ERROR( sourceDistribution, "NULL pointer for source distribution" )
    const Distribution& sourceDist = *sourceDistribution;
    const Communicator& communicator = sourceDist.getCommunicator();
    const IndexType ownID = communicator.getRank();
    const IndexType numPEs = communicator.getSize();
    const IndexType globalN = sourceDist.getGlobalSize();

    mSourceSize = sourceDist.getLocalSize();
    ReadAccess<IndexType> rNewOwners(newOwners);
    SCAI_ASSERT_EQ_ERROR(mSourceSize, rNewOwners.size(), "Got array of " << rNewOwners.size() << " new owners but " << mSourceSize << " local values in distribution.");
    mLocalSourceIndexes = HArray<IndexType>(mSourceSize);

    WriteAccess<IndexType> localSourceIndexes( mLocalSourceIndexes, mSourceSize );
    mNumLocalValues = 0; // count number of local copies from source to target

    for (IndexType i = 0; i < mSourceSize; i++)
    {
        SCAI_ASSERT_DEBUG(rNewOwners[i] < numPEs, "Illegal future owner.");

        if (rNewOwners[i] == ownID)
        {
            localSourceIndexes[mNumLocalValues] = i;
            mNumLocalValues++;
        }
    }

    SCAI_LOG_DEBUG( logger, "Of the " << mSourceSize << " local values, " << mNumLocalValues << " will stay local.")

    //even if numLocalValues == mSourceSize, we cannot know if the distributions are equal, it may have changed for someone else
    localSourceIndexes.resize( mNumLocalValues );

    HaloBuilder::buildFromTargets( newOwners, sourceDist, mHalo );
    const Halo& halo = mHalo;

    //build new list of local indices
    const IndexType haloSize = halo.getHaloSize();
    mTargetSize = halo.getHaloSize() + mNumLocalValues;
    SCAI_ASSERT(communicator.sum(mTargetSize) == globalN, "Sum of target sizes is only " << communicator.sum(mTargetSize) << " of " << globalN)
    SCAI_LOG_DEBUG( logger, "Will get " << halo.getHaloSize() << " values from other PEs for " << mTargetSize << " total.")
    HArray<IndexType> myTargetGlobalIndices(mTargetSize);
    {
        WriteAccess<IndexType> wTargetGlobalIndices(myTargetGlobalIndices);

        for (IndexType i = 0; i < mNumLocalValues; i++)  //maybe improve this with scatter and gather?
        {
            wTargetGlobalIndices[i] = sourceDist.local2global(localSourceIndexes[i]);
        }

        ReadAccess<IndexType> haloIndices(halo.getRequiredIndexes());
        SCAI_ASSERT_EQ_ERROR(haloIndices.size(), haloSize, "Halo inconsistent.");

        for (IndexType i = 0; i < haloSize; i++)
        {
            SCAI_ASSERT_DEBUG(haloIndices[i] < globalN, "Illegal halo index " << haloIndices[i]);
            wTargetGlobalIndices[i+mNumLocalValues]  = haloIndices[i];
        }

        std::sort(wTargetGlobalIndices.get(), wTargetGlobalIndices.get()+mTargetSize);
    }

    mTargetDistribution = DistributionPtr(new GeneralDistribution(globalN, myTargetGlobalIndices, sourceDist.getCommunicatorPtr()));
    const Distribution& targetDist = *mTargetDistribution;
    ContextPtr contextPtr = Context::getHostPtr();
    {
        //assign local target indices
        WriteAccess<IndexType> localTargetIndexes( mLocalTargetIndexes, contextPtr );
        localTargetIndexes.resize(mNumLocalValues);

        for (IndexType i = 0; i < mNumLocalValues; i++)
        {
            IndexType globalI = sourceDist.local2global(localSourceIndexes[i]);
            IndexType targetLocal = targetDist.global2local(globalI);//TODO: maybe optimize
            SCAI_ASSERT_DEBUG( targetLocal < mTargetSize, "Index " << targetLocal << " illegal." );
            localTargetIndexes[i] = targetLocal;
        }
    }

    WriteAccess<IndexType> haloSourceIndexes( mHaloSourceIndexes, contextPtr );
    WriteAccess<IndexType> haloTargetIndexes( mHaloTargetIndexes, contextPtr );
    const CommunicationPlan& providesPlan = halo.getProvidesPlan();
    const CommunicationPlan& requiresPlan = halo.getRequiredPlan();

    SCAI_ASSERT_ERROR(providesPlan.totalQuantity() <= mSourceSize, "Cannot send more indices than I have.");

    ReadAccess<IndexType> haloProvidesIndexes( halo.getProvidesIndexes(), contextPtr );
    ReadAccess<IndexType> haloRequiresIndexes( halo.getRequiredIndexes(), contextPtr );

    //now significant amount of duplicate code. TODO: maybe split off in separate method?
    haloSourceIndexes.resize( providesPlan.totalQuantity() );
    IndexType offset = 0; // runs through halo source indexes

    for ( PartitionId i = 0; i < providesPlan.size(); i++ )
    {
        const IndexType n = providesPlan[i].quantity;
        const IndexType planOffset = providesPlan[i].offset;

        for ( IndexType j = 0; j < n; j++ )
        {
            SCAI_ASSERT_DEBUG( planOffset+j < haloProvidesIndexes.size(), "Index " << planOffset+j << " illegal." );
            haloSourceIndexes[offset++] = haloProvidesIndexes[planOffset+j];
        }
    }

    SCAI_LOG_INFO( logger, "have set " << offset << " halo source indexes" )
    // In contrary to Halo schedules we have here the situation that each non-local
    // index of source should be required by some other processor.
    SCAI_ASSERT_EQ_ERROR( offset, haloSourceIndexes.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mNumLocalValues + offset, mSourceSize, "serious mismatch" )
    // Now add the indexes where to scatter the halo into destination
    haloTargetIndexes.resize( haloSize );
    IndexType targetOffset = 0;

    for ( PartitionId i = 0; i < requiresPlan.size(); i++ )
    {
        const IndexType n = requiresPlan[i].quantity;
        const IndexType planOffset = requiresPlan[i].offset;

        for ( IndexType j = 0; j < n; j++ )
        {
            SCAI_ASSERT_DEBUG( planOffset+j < haloRequiresIndexes.size(), "Index " << planOffset+j << " illegal." );

            haloTargetIndexes[targetOffset] = mTargetDistribution->global2local(haloRequiresIndexes[planOffset+j]);
            SCAI_ASSERT_VALID_INDEX(haloTargetIndexes[targetOffset], mTargetSize, "invalid index");
            targetOffset++;
        }
    }

    SCAI_ASSERT_EQ_ERROR( targetOffset, haloTargetIndexes.size(), "serious mismatch" )

    // all info needed for redistribution is now available
    SCAI_LOG_INFO( logger, "constructed " << *this )
}

/* -------------------------------------------------------------------------- */

void Redistributor::writeAt( std::ostream& stream ) const
{
    stream << "Redistributor( ";
    stream << *mSourceDistribution << "->" << *mTargetDistribution;
    stream << ", " << mSourceSize << "->" << mTargetSize;
    stream << ", local:" << mNumLocalValues;
    stream << ", source halo :" << getHaloSourceSize();
    stream << ", target halo :" << getHaloTargetSize();
    stream << ")";
}

/* -------------------------------------------------------------------------- */

DistributionPtr Redistributor::getSourceDistributionPtr() const
{
    return mSourceDistribution;
}

/* -------------------------------------------------------------------------- */

DistributionPtr Redistributor::getTargetDistributionPtr() const
{
    return mTargetDistribution;
}

/* -------------------------------------------------------------------------- */

void Redistributor::buildVPlans( const IndexType haloSourceSizes[], const IndexType haloTargetSizes[] ) const
{
    const IndexType numProvides = mHalo.getProvidesPlan().totalQuantity();
    const IndexType numRequired = mHalo.getRequiredPlan().totalQuantity();
    // calculate number of provided and required values by summing up the corresponding quantities
    scoped_array<IndexType> provideQuantities( new IndexType[numProvides] );
    scoped_array<IndexType> requiredQuantities( new IndexType[numRequired] );

    // For building the new schedule we need the sizes, can be calculated by the offsets

    for ( IndexType i = 0; i < numProvides; i++ )
    {
        IndexType size = haloSourceSizes[i];
        provideQuantities[i] = size;
        SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
    }

    for ( IndexType i = 0; i < numRequired; i++ )
    {
        IndexType size = haloTargetSizes[i];
        requiredQuantities[i] = size;
        SCAI_LOG_DEBUG( logger, "required[" << i << "] = " << size )
    }

    mProvidesPlan.reset( new CommunicationPlan( mHalo.getProvidesPlan(), provideQuantities.get() ) );
    mRequiredPlan.reset( new CommunicationPlan( mHalo.getRequiredPlan(), requiredQuantities.get() ) );
    SCAI_LOG_INFO( logger, "providesPlan = " << *mProvidesPlan )
    SCAI_LOG_INFO( logger, "requiredPlan = " << *mRequiredPlan )
}

/* -------------------------------------------------------------------------- */

void Redistributor::buildRowPlans(
    const HArray<IndexType>& targetSizes,
    const HArray<IndexType>& sourceSizes ) const
{
    const IndexType numProvides = mHalo.getProvidesPlan().totalQuantity();
    const IndexType numRequired = mHalo.getRequiredPlan().totalQuantity();
    // calculate number of provided and required values by summing up the corresponding quantities
    scoped_array<IndexType> provideQuantities( new IndexType[numProvides] );
    scoped_array<IndexType> requiredQuantities( new IndexType[numRequired] );
    ContextPtr contextPtr = Context::getHostPtr();
    // For building the new schedule we need the sizes, can be calculated by the offsets
    {
        ReadAccess<IndexType> indexes( mHaloSourceIndexes, contextPtr );
        ReadAccess<IndexType> sizes( sourceSizes, contextPtr );

        for ( IndexType i = 0; i < numProvides; i++ )
        {
            IndexType size = sizes[indexes[i]];
            provideQuantities[i] = size;
            SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
        }
    }
    {
        ReadAccess<IndexType> indexes( mHaloTargetIndexes, contextPtr );
        ReadAccess<IndexType> sizes( targetSizes, contextPtr );

        for ( IndexType i = 0; i < numRequired; i++ )
        {
            IndexType size = sizes[indexes[i]];
            requiredQuantities[i] = size;
            SCAI_LOG_DEBUG( logger, "required[" << i << "] = " << size )
        }
    }
    mProvidesPlan.reset( new CommunicationPlan( mHalo.getProvidesPlan(), provideQuantities.get() ) );
    mRequiredPlan.reset( new CommunicationPlan( mHalo.getRequiredPlan(), requiredQuantities.get() ) );
    SCAI_LOG_INFO( logger, "providesPlan = " << *mProvidesPlan )
    SCAI_LOG_INFO( logger, "requiredPlan = " << *mRequiredPlan )
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

/*
 template COMMON_DLL_IMPORTEXPORT
 void Redistributor::redistributeN ( HArray<float>& targetArray,
 const HArray<float>& sourceArray,
 IndexType n ) const;

 template COMMON_DLL_IMPORTEXPORT
 void Redistributor::redistributeN ( HArray<double>& targetArray,
 const HArray<double>& sourceArray,
 IndexType n ) const;
 */

} /* end namespace dmemo */

} /* end namespace scai */
