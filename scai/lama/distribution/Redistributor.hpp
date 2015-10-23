/**
 * @file Redistributor.hpp
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
 * @brief Implementation of a class that handles redistribution of matrices and vectors
 * @author Thomas Brandes
 * @date 08.10.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/lama/distribution/Distribution.hpp>
#include <scai/lama/distribution/Halo.hpp>
#include <scai/lama/LAMAKernel.hpp>
#include <scai/lama/UtilsInterface.hpp>

// internal scai libraries
#include <scai/hmemo/LAMAArray.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/unique_ptr.hpp>

namespace scai
{

using namespace scai::tasking;

namespace lama
{

using hmemo::LAMAArray;

/** This class allows to create objects that handle redistributions of vector and
 *  matrices from one distribution into another distribution.
 *
 *  Once created, it has built internal data structures like communicaton plans
 *  that restrict the redistribution just to the transfer of the corresponding data.
 */

class COMMON_DLL_IMPORTEXPORT Redistributor: public scai::common::Printable
{
public:

    /** Build an object by the source and target distribution.
     *
     *  @param[in] targetDistribution  the new distribution
     *  @param[in] sourceDistribution  the old distribution
     *
     *  The global size of both distributions must be the same.
     */

    Redistributor( DistributionPtr targetDistribution, DistributionPtr sourceDistribution );

    /** Getter needed for distributions */

    DistributionPtr getTargetDistributionPtr() const;

    DistributionPtr getSourceDistributionPtr() const;

    IndexType getSourceLocalSize() const
    {
        return mSourceSize;
    }

    IndexType getTargetLocalSize() const
    {
        return mTargetSize;
    }

    /** Redistribution of a distributed vector as LAMAArrays.
     *
     *  @param[out] targetArray  vector in target distribution
     *  @param[in]  sourceArray  vector in source distribution
     *
     *  Size of sourceArray must be at least sourceDistribution->getLocalSize()
     *  Size of targetArray must be at least targetDistribution->getLocalSize()
     *
     *  This routine is a simple straightforward solution. More
     *  efficient implementations are possible: overlapping local copies  and
     *  communication.
     */
    template<typename ValueType>
    void redistribute( LAMAArray<ValueType>& targetArray, const LAMAArray<ValueType>& sourceArray ) const;

    /** Redistribution of a distributed vector as LAMAArrays.
     *
     *  @param[out] targetArray  vector in target distribution
     *  @param[in]  sourceArray  vector in source distribution
     *  @param[in]  n            TODO[doxy] Complete Description.
     *
     *  Size of sourceArray must be at least sourceDistribution->getLocalSize()
     *  Size of targetArray must be at least targetDistribution->getLocalSize()
     *
     *  This code only demonstrates how the redistributor schedule works. More
     *  efficient implementations are possible: overlapping local copies  and
     *  communication.
     */
    template<typename ValueType>
    void redistributeN( LAMAArray<ValueType>& targetArray, const LAMAArray<ValueType>& sourceArray, IndexType n ) const;

    /** Redistribution of ragged arrays. */

    template<typename ValueType>
    void redistributeV(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetOffsets,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceOffsets ) const;

    template<typename ValueType>
    static void gather(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceIndexes )
    {
        using namespace scai::hmemo;

        static LAMAKernel<UtilsInterface::setGather<ValueType, ValueType> > setGather;

        ContextPtr loc = Context::getHostPtr();   // do it on host

        IndexType n = sourceIndexes.size();

        WriteOnlyAccess<ValueType> target( targetArray, loc, n );
        ReadAccess<ValueType> source( sourceArray, loc );
        ReadAccess<IndexType> indexes( sourceIndexes, loc );

        setGather[loc]( target.get(), source.get(), indexes.get(), indexes.size() );
    }

    template<typename ValueType>
    static void gatherN(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceIndexes,
        const IndexType n )
    {
        using namespace scai::hmemo;

        ContextPtr loc = Context::getHostPtr();

        WriteAccess<ValueType> target( targetArray, loc );
        ReadAccess<ValueType> source( sourceArray, loc );
        ReadAccess<IndexType> indexes( sourceIndexes, loc );

        #pragma omp parallel for

        for( IndexType i = 0; i < indexes.size(); i++ )
        {
            SCAI_LOG_DEBUG( logger,
                            "targetN[" << i << "] = sourceN[" << indexes[i] << "] = " << source[indexes[i] * n] << " ..." )

            for( IndexType j = 0; j < n; j++ )
            {
                target[i * n + j] = source[indexes[i] * n + j];
            }
        }
    }

    template<typename ValueType>
    static void gatherV(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceOffsets,
        const LAMAArray<IndexType>& sourceIndexes );

    template<typename ValueType>
    static void scatter(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray )
    {
        using namespace scai::hmemo;

        ContextPtr loc = Context::getHostPtr();

        WriteAccess<ValueType> target( targetArray, loc );
        ReadAccess<IndexType> indexes( targetIndexes, loc );
        ReadAccess<ValueType> source( sourceArray, loc );

        for( IndexType i = 0; i < indexes.size(); i++ )
        {
            SCAI_LOG_DEBUG( logger, "target[" << indexes[i] << "] = source[" << i << "] = " << source[i] )

            target[indexes[i]] = source[i];
        }
    }

    template<typename ValueType>
    static void scatterN(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray,
        const IndexType n )
    {
        using namespace scai::hmemo;

        ContextPtr loc = Context::getHostPtr();

        WriteAccess<ValueType> target( targetArray, loc );
        ReadAccess<IndexType> indexes( targetIndexes, loc );
        ReadAccess<ValueType> source( sourceArray, loc );

        #pragma omp parallel for

        for( IndexType i = 0; i < indexes.size(); i++ )
        {
            SCAI_LOG_DEBUG( logger,
                            "targetN[" << indexes[i] << "] = sourceN[" << i << "] = " << source[i * n] << " ..." )

            for( IndexType j = 0; j < n; j++ )
            {
                target[indexes[i] * n + j] = source[i * n + j];
            }
        }
    }

    template<typename ValueType>
    static void scatterV(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetOffsets,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray );

    template<typename ValueType>
    static void copy(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceIndexes )
    {
        using namespace scai::hmemo;

        ContextPtr loc = Context::getHostPtr();

        WriteAccess<ValueType> target( targetArray, loc );
        ReadAccess<ValueType> source( sourceArray, loc );
        ReadAccess<IndexType> tindexes( targetIndexes, loc );
        ReadAccess<IndexType> sindexes( sourceIndexes, loc );

        SCAI_ASSERT_ERROR( tindexes.size() == sindexes.size(), "index size mismatch" )

        for( IndexType i = 0; i < tindexes.size(); i++ )
        {
            SCAI_LOG_DEBUG( logger,
                            "target[" << tindexes[i] << "] = source[" << sindexes[i] << "] = " << source[ sindexes[i] ] )

            target[tindexes[i]] = source[sindexes[i]];
        }
    }

    template<typename ValueType>
    static void copyN(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceIndexes,
        IndexType n )
    {
        using namespace scai::hmemo;

        ContextPtr loc = Context::getHostPtr();

        WriteAccess<ValueType> target( targetArray, loc );
        ReadAccess<ValueType> source( sourceArray, loc );
        ReadAccess<IndexType> tindexes( targetIndexes, loc );
        ReadAccess<IndexType> sindexes( sourceIndexes, loc );

        SCAI_ASSERT_ERROR( tindexes.size() == sindexes.size(), "index size mismatch" )

        #pragma omp parallel for

        for( IndexType i = 0; i < tindexes.size(); i++ )
        {
            SCAI_LOG_DEBUG( logger,
                            "targetN[" << tindexes[i] << "] = sourceN[" << sindexes[i] << "] = " << source[ sindexes[i] * n ] << " ..." )

            for( IndexType j = 0; j < n; j++ )
            {
                target[tindexes[i] * n + j] = source[sindexes[i] * n + j];
            }
        }
    }

    template<typename ValueType>
    static void copyV(
        LAMAArray<ValueType>& targetArray,
        const LAMAArray<IndexType>& targetOffsets,
        const LAMAArray<IndexType>& targetIndexes,
        const LAMAArray<ValueType>& sourceArray,
        const LAMAArray<IndexType>& sourceOffsets,
        const LAMAArray<IndexType>& sourceIndexes );

    IndexType getHaloSourceSize() const
    {
        return mHaloSourceIndexes.size();
    }
    IndexType getHaloTargetSize() const
    {
        return mHaloTargetIndexes.size();
    }

    template<typename ValueType>
    void exchangeHalo( LAMAArray<ValueType>& targetHalo, const LAMAArray<ValueType>& sourceHalo ) const;

    template<typename ValueType>
    void exchangeHaloN(
        LAMAArray<ValueType>& targetHalo,
        const LAMAArray<ValueType>& sourceHalo,
        const IndexType n ) const;

    void buildVPlans( const IndexType haloSourceSizes[], const IndexType haloTargetSizes[] ) const;

    /** The redistributor can also be used to exchange rows of distributed sparse matrices instead
     *  of distributed vector elements. This method will build the corresponding exchange schedule.
     */

    void buildRowPlans( const LAMAArray<IndexType>& targetSizes, const LAMAArray<IndexType>& sourceSizes ) const;

    IndexType getVHaloSourceSize() const
    {
        return mProvidesPlan->totalQuantity();
    }
    IndexType getVHaloTargetSize() const
    {
        return mRequiredPlan->totalQuantity();
    }

    template<typename ValueType>
    void exchangeVHalo( LAMAArray<ValueType>& targetHalo, const LAMAArray<ValueType>& sourceHalo ) const;

    const LAMAArray<IndexType>& getLocalSourceIndexes() const
    {
        return mLocalSourceIndexes;
    }
    ;
    const LAMAArray<IndexType>& getLocalTargetIndexes() const
    {
        return mLocalTargetIndexes;
    }
    ;
    const LAMAArray<IndexType>& getHaloSourceIndexes() const
    {
        return mHaloSourceIndexes;
    }
    ;
    const LAMAArray<IndexType>& getHaloTargetIndexes() const
    {
        return mHaloTargetIndexes;
    }
    ;

private:

    virtual void writeAt( std::ostream& stream ) const;

    DistributionPtr mSourceDistribution;
    DistributionPtr mTargetDistribution;

    IndexType mSourceSize; // = mSourceDistribution->getLocalSize()
    IndexType mTargetSize; // = mTargetDistribution->getLocalSize()

    LAMAArray<IndexType> mLocalSourceIndexes;
    LAMAArray<IndexType> mLocalTargetIndexes;

    LAMAArray<IndexType> mHaloSourceIndexes;
    LAMAArray<IndexType> mHaloTargetIndexes;

    IndexType mNumLocalValues; // common number of local values

    Halo mHalo; // Halo structure for exchanging non-local values

    mutable common::unique_ptr<CommunicationPlan> mProvidesPlan;
    mutable common::unique_ptr<CommunicationPlan> mRequiredPlan;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::redistribute( LAMAArray<ValueType>& targetArray, const LAMAArray<ValueType>& sourceArray ) const
{
    using namespace scai::hmemo;

    SCAI_REGION( "Redistributor.redistribute" )

    {
        // make sure that target array has sufficient memory

        WriteOnlyAccess<ValueType> target( targetArray, mTargetSize );
    }

    // allocate memory for source (provides) and target (required) halo

    LAMAArray<ValueType> sourceHalo( getHaloSourceSize() );
    LAMAArray<ValueType> targetHalo( getHaloTargetSize() );

    SCAI_LOG_DEBUG( logger, "gather: sourceHalo " << mHaloSourceIndexes.size() << " values" )

    gather( sourceHalo, sourceArray, mHaloSourceIndexes );

    SCAI_LOG_DEBUG( logger, "copy: source -> target " << mLocalTargetIndexes.size() << " values" )

    copy( targetArray, mLocalTargetIndexes, sourceArray, mLocalSourceIndexes );

    exchangeHalo( targetHalo, sourceHalo );

    SCAI_LOG_DEBUG( logger, "scatter: targetHalo " << mHaloTargetIndexes.size() << " values" )

    scatter( targetArray, mHaloTargetIndexes, targetHalo );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::redistributeN(
    LAMAArray<ValueType>& targetArray,
    const LAMAArray<ValueType>& sourceArray,
    IndexType n ) const
{
    using namespace scai::hmemo;

    SCAI_REGION( "Redistributor.redistributeN" )

    ContextPtr loc = Context::getHostPtr();

    {
        // make sure that target array has sufficient memory

        WriteOnlyAccess<ValueType> target( targetArray, loc, mTargetSize * n );
    }

    // allocate memory for source (provides) and target (required) halo

    LAMAArray<ValueType> sourceHalo( n * getHaloSourceSize() );
    LAMAArray<ValueType> targetHalo( n * getHaloTargetSize() );

    SCAI_LOG_DEBUG( logger, "gather: sourceHalo " << mHaloSourceIndexes.size() << " * " << n << " values" )

    gatherN( sourceHalo, sourceArray, mHaloSourceIndexes, n );

    SCAI_LOG_DEBUG( logger, "copy: source -> target " << mLocalTargetIndexes.size() << " * " << n << " values" )

    copyN( targetArray, mLocalTargetIndexes, sourceArray, mLocalSourceIndexes, n );

    exchangeHaloN( targetHalo, sourceHalo, n );

    SCAI_LOG_DEBUG( logger, "scatter: targetHalo " << mHaloTargetIndexes.size() << " * " << n << " values" )

    scatterN( targetArray, mHaloTargetIndexes, targetHalo, n );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::redistributeV(
    LAMAArray<ValueType>& targetArray,
    const LAMAArray<IndexType>& targetOffsets,
    const LAMAArray<ValueType>& sourceArray,
    const LAMAArray<IndexType>& sourceOffsets ) const
{
    SCAI_REGION( "Redistributor.redistributeV" )

    // allocate memory for source (provides) and target (required) halo

    LAMAArray<ValueType> sourceHalo( getVHaloSourceSize() );
    LAMAArray<ValueType> targetHalo( getVHaloTargetSize() );

    gatherV( sourceHalo, sourceArray, sourceOffsets, getHaloSourceIndexes() );

    copyV( targetArray, targetOffsets, mLocalTargetIndexes, sourceArray, sourceOffsets, mLocalSourceIndexes );

    exchangeVHalo( targetHalo, sourceHalo );

    scatterV( targetArray, targetOffsets, mHaloTargetIndexes, targetHalo );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::gatherV(
    LAMAArray<ValueType>& targetArray,
    const LAMAArray<ValueType>& sourceArray,
    const LAMAArray<IndexType>& sourceOffsets,
    const LAMAArray<IndexType>& sourceIndexes )
{
    using namespace scai::hmemo;

    const IndexType n = sourceIndexes.size();

    ContextPtr loc = Context::getHostPtr();

    WriteAccess<ValueType> wTargetArray( targetArray, loc );
    ReadAccess<ValueType> rSourceArray( sourceArray, loc );
    ReadAccess<IndexType> rSourceOffsets( sourceOffsets, loc );
    ReadAccess<IndexType> rSourceIndexes( sourceIndexes, loc );

    // Note: we have no target offsets array

    IndexType targetOffset = 0;

    for( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = rSourceIndexes[ii];

        for( IndexType j = rSourceOffsets[i]; j < rSourceOffsets[i + 1]; ++j )
        {
            wTargetArray[targetOffset++] = rSourceArray[j];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::scatterV(
    LAMAArray<ValueType>& targetArray,
    const LAMAArray<IndexType>& targetOffsets,
    const LAMAArray<IndexType>& targetIndexes,
    const LAMAArray<ValueType>& sourceArray )
{
    using namespace scai::hmemo;

    ContextPtr loc = Context::getHostPtr();

    const IndexType n = targetIndexes.size();

    WriteAccess<ValueType> wTargetArray( targetArray, loc );
    ReadAccess<IndexType> rTargetOffsets( targetOffsets, loc );
    ReadAccess<IndexType> rTargetIndexes( targetIndexes, loc );
    ReadAccess<ValueType> rSourceArray( sourceArray, loc );

    // Note: we have no source offsets array, no parallelization possible

    IndexType sourceOffset = 0;

    for( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = rTargetIndexes[ii];

        for( IndexType j = rTargetOffsets[i]; j < rTargetOffsets[i + 1]; ++j )
        {
            wTargetArray[j] = rSourceArray[sourceOffset++];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::copyV(
    LAMAArray<ValueType>& targetArray,
    const LAMAArray<IndexType>& targetOffsets,
    const LAMAArray<IndexType>& targetIndexes,
    const LAMAArray<ValueType>& sourceArray,
    const LAMAArray<IndexType>& sourceOffsets,
    const LAMAArray<IndexType>& sourceIndexes )
{
    using namespace scai::hmemo;

    SCAI_ASSERT_EQUAL_ERROR( targetIndexes.size(), sourceIndexes.size() )

    ContextPtr loc = Context::getHostPtr();

    const IndexType n = targetIndexes.size();

    WriteAccess<ValueType> wTargetArray( targetArray, loc );
    ReadAccess<IndexType> rTargetOffsets( targetOffsets, loc );
    ReadAccess<IndexType> rTargetIndexes( targetIndexes, loc );
    ReadAccess<ValueType> rSourceArray( sourceArray, loc );
    ReadAccess<IndexType> rSourceOffsets( sourceOffsets, loc );
    ReadAccess<IndexType> rSourceIndexes( sourceIndexes, loc );

    for( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType sourceI = rSourceIndexes[ii];
        IndexType targetI = rTargetIndexes[ii];

        IndexType k = rTargetOffsets[targetI];

        for( IndexType j = rSourceOffsets[sourceI]; j < rSourceOffsets[sourceI + 1]; ++j )
        {
            wTargetArray[k] = rSourceArray[j];
            ++k;
        }

        SCAI_ASSERT_EQUAL_DEBUG( k, rTargetOffsets[targetI + 1] )
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::exchangeHalo( LAMAArray<ValueType>& targetHalo, const LAMAArray<ValueType>& sourceHalo ) const
{
    SCAI_REGION( "Redistributor.exchangeHalo" )

    const Communicator& comm = mSourceDistribution->getCommunicator();

    // use asynchronous communication to avoid deadlocks

    SyncToken* token = comm.exchangeByPlanAsync( targetHalo, mHalo.getRequiredPlan(), sourceHalo,
                       mHalo.getProvidesPlan() );

    token->wait();

    delete token;

    // synchronization is done implicitly
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::exchangeHaloN(
    LAMAArray<ValueType>& targetHalo,
    const LAMAArray<ValueType>& sourceHalo,
    const IndexType n ) const
{
    SCAI_REGION( "Redistributor.exchangeHaloN" )

    const Communicator& comm = mSourceDistribution->getCommunicator();

    // Communication plans are built by multiplication with n

    CommunicationPlan requiredN( mHalo.getRequiredPlan(), n );
    CommunicationPlan providesN( mHalo.getProvidesPlan(), n );

    SCAI_LOG_DEBUG( logger, "requiredN ( n = " << n << "): " << requiredN )
    SCAI_LOG_DEBUG( logger, "providesN ( n = " << n << "): " << providesN )

    // use asynchronous communication to avoid deadlocks

    comm.exchangeByPlan( targetHalo, requiredN, sourceHalo, providesN );

    // synchronization is done implicitly at the end of this scope
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void Redistributor::exchangeVHalo( LAMAArray<ValueType>& targetHalo, const LAMAArray<ValueType>& sourceHalo ) const
{
    SCAI_REGION( "Redistributor.exchangeVHalo" )

    const Communicator& comm = mSourceDistribution->getCommunicator();

    SCAI_ASSERT_ERROR( mRequiredPlan.get(), "There was no previous call of buildVPlan" )

    delete comm.exchangeByPlanAsync( targetHalo, *mRequiredPlan, sourceHalo, *mProvidesPlan );
}

} /* end namespace lama */

} /* end namespace scai */
