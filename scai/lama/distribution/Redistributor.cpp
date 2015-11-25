/**
 * @file Redistributor.cpp
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
 * @brief Implementation of methods for Redistributor class.
 * @author Thomas Brandes
 * @date 08.10.2011
 */

// hpp
#include <scai/lama/distribution/Redistributor.hpp>

// local library
#include <scai/lama/distribution/HaloBuilder.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::unique_ptr;
using common::scoped_array;

namespace lama
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

    SCAI_ASSERT_EQUAL_ERROR( sourceDist.getGlobalSize(), targetDist.getGlobalSize() )
    SCAI_ASSERT_EQUAL_ERROR( sourceDist.getCommunicator(), targetDist.getCommunicator() )

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

    for( IndexType i = 0; i < mTargetSize; i++ )
    {
        IndexType globalIndex = targetDist.local2global( i );

        if( sourceDist.isLocal( globalIndex ) )
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

    HaloBuilder::build( sourceDist, requiredIndexes, mHalo );

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

    for( IndexType i = 0; i < providesPlan.size(); i++ )
    {
        IndexType n = providesPlan[i].quantity;

        const IndexType* pindexes = haloProvidesIndexes.get() + providesPlan[i].offset;

        for( IndexType j = 0; j < n; j++ )
        {
            SCAI_LOG_TRACE( logger, "halo source index[" << offset << "] = " << pindexes[j] )

            haloSourceIndexes[offset++] = pindexes[j];
        }
    }

    SCAI_LOG_INFO( logger, "have set " << offset << " halo source indexes" )

    // In contrary to Halo schedules we have here the situation that each non-local
    // index of source should be required by some other processor.

    SCAI_ASSERT_EQUAL_ERROR( offset, haloSourceIndexes.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumLocalValues + offset, mSourceSize )

    // Now add the indexes where to scatter the halo into destination

    IndexType haloSize = halo.getHaloSize();

    SCAI_ASSERT_ERROR( mNumLocalValues + haloSize == mTargetSize, "size mismatch" )

    haloTargetIndexes.resize( haloSize );

    for( IndexType i = 0; i < haloSize; i++ )
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

    for( IndexType i = 0; i < numProvides; i++ )
    {
        IndexType size = haloSourceSizes[i];
        provideQuantities[i] = size;
        SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
    }

    for( IndexType i = 0; i < numRequired; i++ )
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

        for( IndexType i = 0; i < numProvides; i++ )
        {
            IndexType size = sizes[indexes[i]];
            provideQuantities[i] = size;
            SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
        }
    }

    {
        ReadAccess<IndexType> indexes( mHaloTargetIndexes, contextPtr );
        ReadAccess<IndexType> sizes( targetSizes, contextPtr );

        for( IndexType i = 0; i < numRequired; i++ )
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

} /* end namespace lama */

} /* end namespace scai */
