/**
 * @file GPICommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief GPICommunicator.hpp
 * @author Lauretta Schubert, Thomas Brandes
 * @date 25.02.2014
 * @since 1.1.0
 */

#pragma once

#include <GASPI.h>

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/CRTPCommunicator.hpp>
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

class COMMON_DLL_IMPORTEXPORT GPICommunicator: public CRTPCommunicator<GPICommunicator>, 
                                               public Communicator::Register<GPICommunicator>,
                                               public common::enable_shared_from_this<GPICommunicator>

{
    // Only GPICommunicatorManager is allowed to create GPI communicator

    friend class GPISyncToken;
    friend class CRTPCommunicator<GPICommunicator>;

public:

    virtual ~GPICommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    /** @brief Provide implementation for Communicator::getSize */

    virtual PartitionId getSize() const;

    /** @brief Provide implementation for Communicator::getRank */

    virtual PartitionId getRank() const;

    /** @brief Provide implementation for Communicator::getNodeSize */

    virtual PartitionId getNodeSize() const;

    /** @brief Provide implementation for Communicator::getNodeRank */

    virtual PartitionId getNodeRank() const;

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

    void setNodeData();

    template<typename T>
    inline static gaspi_datatype_t getGPIType();

    template<typename T>
    T sumImpl( T myVal ) const;
    template<typename T>
    T maxImpl( T myVal ) const;
    template<typename T>
    T minImpl( T myVal ) const;

    template<typename T>
    void bcastImpl( T val[], const IndexType n, const PartitionId root ) const;

    template<typename T>
    void scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const;

    template<typename T>
    void scatterVImpl(
        T myvals[],
        const IndexType n,
        const PartitionId root,
        const T allvals[],
        const IndexType sizes[] ) const;

    template<typename T>
    void gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const;

    template<typename T>
    void gatherVImpl(
        T allvals[],
        const IndexType n,
        const PartitionId root,
        const T myvals[],
        const IndexType sizes[] ) const;

    template<typename T>
    IndexType shiftImpl(
        T newvals[],
        const IndexType newSize,
        const PartitionId source,
        const T oldVals[],
        const IndexType oldSize,
        const PartitionId dest ) const;

    template<typename T>
    tasking::SyncToken* shiftAsyncImpl(
        T newvals[],
        const PartitionId source,
        const T oldVals[],
        const PartitionId dest,
        const IndexType size ) const;

    template<typename T>
    void swapImpl( T val[], const IndexType n, PartitionId partner ) const;

    template<typename T>
    void maxlocImpl( T& val, IndexType& location, PartitionId root ) const;

    template<typename T>
    void exchangeByPlanImpl(
        T recvData[],
        const CommunicationPlan& recvPlan,
        const T sendData[],
        const CommunicationPlan& sendPlan ) const;

    template<typename T>
    tasking::SyncToken* exchangeByPlanAsyncImpl(
        T recvData[],
        const CommunicationPlan& recvPlan,
        const T sendData[],
        const CommunicationPlan& sendPlan ) const;

    template<typename ValueType>
    void all2allvImpl( ValueType* recvBuffer[], IndexType recvCount[], ValueType* sendBuffer[], IndexType sendCount[] ) const;

    const gaspi_queue_id_t mQueueID;

    gaspi_rank_t mRank; // rank of this processor
    gaspi_rank_t mSize; // size of communicator

    // Be careful here: remoteOffset must be set to ptr of remSegment on other processor */

    template<typename T>
    void remoteWrite( const SegmentData<T>& localSegment, const IndexType localOffset, const PartitionId remP,
                      SegmentData<T>& remSegment, const IndexType remOffset, const IndexType size ) const;

    // Be careful here: dstOffset must be the logical offset

    template<typename T>
    void localWrite( const SegmentData<T>& srcSegment, const IndexType srcOffset, 
                     SegmentData<T>& dstSegment, const IndexType dstOffset, const IndexType size ) const;

    template<typename T>
    void remoteRead( SegmentData<T>& localSegment, const IndexType localOffset, const PartitionId remP,
                     const SegmentData<T>& remSegment, const IndexType remOffset, const IndexType size ) const;

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

