/**
 * @file GPICommunicator.hpp
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
 * @brief GPICommunicator.hpp
 * @author Lauretta Schubert, Thomas Brandes
 * @date 25.02.2014
 */

#pragma once

#include <GASPI.h>

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/gpi/SegmentData.hpp>

// others
#include <scai/tasking/SyncToken.hpp>

// logging
#include <scai/logging.hpp>

// common
#include <scai/common/shared_ptr.hpp>
#include <scai/common/SCAITypes.hpp>

#include <vector>

namespace scai
{

namespace dmemo
{

/** Communicator class that implements communication and data exchange via GPI.
 *
 *  gaspi_proc_init is called in the constructor, gaspi_proc_term is called in the destructor.
 */

class COMMON_DLL_IMPORTEXPORT GPICommunicator: public Communicator,
    public Communicator::Register<GPICommunicator>,
    public common::enable_shared_from_this<GPICommunicator>

{
    // Only GPICommunicatorManager is allowed to create GPI communicator

    friend class GPISyncToken;

public:

    virtual ~GPICommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    /** All-to-all exchange of an integer value between all processors.
     *
     * @param[out] recvValues will contain one value from each processor
     * @param[in]  sendValues must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( IndexType recvValues[], const IndexType sendValues[] ) const;

    virtual void synchronize() const;

    /** scatter */

    virtual void writeAt( std::ostream& stream ) const;

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory

    static CommunicatorKind createValue();

private:

    GPICommunicator( );

    void wait() const;

    /** Implementation of Communicator::getProcessorName */

    virtual void getProcessorName( char* name ) const;

    /** Implementation of Communicator::maxProcessorName */

    virtual size_t maxProcessorName() const;

    inline static gaspi_datatype_t getGPIType( common::scalar::ScalarType stype );

    /** Implementation of Communicator::sumImpl */

    virtual void sumImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::minImpl */

    virtual void minImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::maxImpl */

    virtual void maxImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation for pure method Communicator::bcastImpl */

    void bcastImpl( void* val, const IndexType n, const PartitionId root, common::scalar::ScalarType stype ) const;


    /** Implementation of pure method Communicator::scatterImpl */

    void scatterImpl( void* myVals, const IndexType n, const PartitionId root, const void* allVals, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::scatterVImpl */

    void scatterVImpl( void* myVals, const IndexType n, const PartitionId root,
                       const void* allVals, const IndexType sizes[], common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherImpl */

    void gatherImpl( void* allVals, const IndexType n, const PartitionId root, const void* myVals, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherVImpl */

    void gatherVImpl(
        void* allvals,
        const IndexType n,
        const PartitionId root,
        const void* myvals,
        const IndexType sizes[],
        common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::shiftImpl */

    IndexType shiftImpl(
        void* newVals,
        const IndexType newSize,
        const PartitionId source,
        const void* oldVals,
        const IndexType oldSize,
        const PartitionId dest,
        common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::shiftAsyncImpl */

    tasking::SyncToken* shiftAsyncImpl(
        void* newVals,
        const PartitionId source,
        const void* oldVals,
        const PartitionId dest,
        const IndexType size,
        common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::swapImpl */

    void swapImpl( void* val, const IndexType n, PartitionId partner, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::supportsLocReduction */

    virtual bool supportsLocReduction( common::scalar::ScalarType vType, common::scalar::ScalarType iType ) const;

    void maxlocImpl( void* val, IndexType* location, PartitionId root, common::scalar::ScalarType stype ) const;

    void minlocImpl( void* val, IndexType* location, PartitionId root, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::exchangeByPlanImpl */

    void exchangeByPlanImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::exchangeByPlanAsyncImpl */

    tasking::SyncToken* exchangeByPlanAsyncImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::scalar::ScalarType stype ) const;

    /** GPI Implementation for pure method Communciator::all2allvImpl */

    void all2allvImpl( void* recvBuffer[], IndexType recvCount[],
                       void* sendBuffer[], IndexType sendCount[],
                       common::scalar::ScalarType stype ) const;

    const gaspi_queue_id_t mQueueID;

    gaspi_rank_t gRank; // rank of this processor (same as mRank, but might be other type)
    gaspi_rank_t gSize; // size of communicator

    // Be careful here: remoteOffset must be set to ptr of remSegment on other processor */

    void remoteWrite( const SegmentData& localSegment, const IndexType localOffset, const PartitionId remP,
                      SegmentData& remSegment, const IndexType remOffset, const IndexType size ) const;

    // Be careful here: dstOffset must be the logical offset

    void localWrite( const SegmentData& srcSegment, const IndexType srcOffset,
                     SegmentData& dstSegment, const IndexType dstOffset, const IndexType size ) const;

    void remoteRead( SegmentData& localSegment, const IndexType localOffset, const PartitionId remP,
                     const SegmentData& remSegment, const IndexType remOffset, const IndexType size ) const;

    /** Helper routine regarding notification. */
    void notify( const gaspi_segment_id_t segID,
                 PartitionId target,
                 PartitionId pos,
                 IndexType val ) const;

    IndexType notifyWait( const gaspi_segment_id_t segID, PartitionId pos ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static const gaspi_timeout_t timeout; // = GASPI_BLOCK
    static const gaspi_group_t   group;   // = GASPI_GROUP_ALL

    static const int defaultTag;

    //!< maximal number of notifications, should be >= 2 * mSize
    PartitionId mMaxNotifications;

protected:

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

private:

    /** Common implementation for user reduction */

    void reduce( void* outValues, const void* inValues, const IndexType n, gaspi_reduce_operation_t op, gaspi_size_t elem_size ) const;

    // Guard class whose destructor takes care of GPI finalize

    class GPIGuard
    {
    public:

        GPIGuard();

        ~GPIGuard();
    };

    static GPIGuard guard;   // define one guard variable
};

} // namespace dmemo

} // namespace scai

