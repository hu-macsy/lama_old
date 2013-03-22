/**
 * @file PGASCommunicator.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief PGASCommunicator.hpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */
#ifndef LAMA_PGASCOMMUNICATOR_H_
#define LAMA_PGASCOMMUNICATOR_H_

#include <lama/Communicator.hpp>

#include <lama/pgas/PGASInterface.hpp>

namespace lama
{

class PGASCommunicator: public lama::Communicator
{

    // Only PGASCommunicatorManager is allowed to create PGAS communicator

    friend class PGASCommunicatorManager;
    PGASCommunicator( int& argc, char** & argv );

public:
    virtual ~PGASCommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    virtual PartitionId getSize() const;

    virtual PartitionId getRank() const;

    /** All-to-all exchange of an integer value between all processors.
     *
     * @param[out]   recvValues will contain one value from each processor
     * @param[in]    sendValues must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( int* recvValues, const int* sendValues ) const;

    /** Exchange data between all processors by communication plans.
     *
     *  @param[out] recvData   buffer for data received from other processors
     *  @param[in]  recvPlan   number of elements and offsets for receiving
     *  @param[in]  sendData   buffer for data to send to other processors
     *  @param[in]  sendPlan   contains number of elements and offsets for sending data
     */
    virtual void exchangeByPlan(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual void exchangeByPlan(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual void exchangeByPlan(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual std::auto_ptr<SyncToken> exchangeByPlanAsync(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual std::auto_ptr<SyncToken> exchangeByPlanAsync(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual std::auto_ptr<SyncToken> exchangeByPlanAsync(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual float sum( const float value ) const;

    virtual double sum( const double value ) const;

    virtual int sum( const int value ) const;

    virtual float min( const float value ) const;

    virtual size_t sum( const size_t value ) const;

    virtual float max( const float value ) const;

    virtual double min( const double value ) const;

    virtual double max( const double value ) const;

    virtual int min( const int value ) const;

    virtual int max( const int value ) const;

    virtual void gather( std::vector<float>& values, float value ) const;

    virtual void synchronize() const;

    virtual IndexType shift(
        double newVals[],
        const IndexType newSize,
        const double oldVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual IndexType shift(
        float newVals[],
        const IndexType newSize,
        const float oldVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual IndexType shift(
        int newVals[],
        const IndexType newSize,
        const int oldVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual std::auto_ptr<SyncToken> shiftAsync(
        double newVals[],
        const double oldVals[],
        const IndexType size,
        const int direction ) const;

    virtual std::auto_ptr<SyncToken> shiftAsync(
        float newVals[],
        const float oldVals[],
        const IndexType size,
        const int direction ) const;

    virtual std::auto_ptr<SyncToken> shiftAsync(
        int newVals[],
        const int oldVals[],
        const IndexType size,
        const int direction ) const;

    /** Broadcast an array of doubles from root to all processors in simulation */

    virtual void bcast( double val[], const IndexType n, const PartitionId root ) const;

    /** Broadcast an array of int from root to all processors in simulation */

    virtual void bcast( int val[], const IndexType n, const PartitionId root ) const;

    /** Broadcast an array of float from root to all processors in simulation */

    virtual void bcast( float val[], const IndexType n, const PartitionId root ) const;

    /** scatter */

    virtual void scatter( double myvals[], const IndexType n, const PartitionId root, const double allvals[] ) const;
    virtual void scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const;
    virtual void scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const;

    virtual void scatter(
        double myvals[],
        const IndexType n,
        const PartitionId root,
        const double allvals[],
        const IndexType sizes[] ) const;
    virtual void scatter(
        float myvals[],
        const IndexType n,
        const PartitionId root,
        const float allvals[],
        const IndexType sizes[] ) const;
    virtual void scatter(
        int myvals[],
        const IndexType n,
        const PartitionId root,
        const int allvals[],
        const IndexType sizes[] ) const;

    /** gather */

    virtual void gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const;
    virtual void gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const;
    virtual void gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const;

    virtual void gather(
        double allvals[],
        const IndexType n,
        const PartitionId root,
        const double myvals[],
        const IndexType sizes[] ) const;
    virtual void gather(
        float allvals[],
        const IndexType n,
        const PartitionId root,
        const float myvals[],
        const IndexType sizes[] ) const;
    virtual void gather(
        int allvals[],
        const IndexType n,
        const PartitionId root,
        const int myvals[],
        const IndexType sizes[] ) const;

    /** Maximal value combined with a location value where maximum was found
     *
     *  @param[in,out] val        is a value on each processor, only out for root with maximal value
     *  @param[in,out] location   is an additional int value, only out for root
     *  @param[in]     root       TODO[doxy] Complete Description.
     *
     *  Only root processor will contain the maximal value and the location loc.
     */
    virtual void maxloc( double& val, int& location, const PartitionId root ) const;

    virtual void maxloc( float& val, int& location, const PartitionId root ) const;

    virtual void maxloc( int& val, int& location, const PartitionId root ) const;

    /** This routine swaps an array with a partner processor. */

    virtual void swap( double val[], const IndexType n, const PartitionId partner ) const;

    virtual void swap( float val[], const IndexType n, const PartitionId partner ) const;

    virtual void swap( int val[], const IndexType n, const PartitionId partner ) const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    template<typename T>
    inline void send( const T* buffer, int count, int target ) const;

    template<typename T>
    void scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const;

    template<typename T>
    void scatterImpl(
        T myvals[],
        const IndexType n,
        const PartitionId root,
        const T allvals[],
        const IndexType sizes[] ) const;

    template<typename T>
    void gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const;

    template<typename T>
    void gatherImpl(
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
    std::auto_ptr<SyncToken> shiftAsyncImpl(
        T newvals[],
        const PartitionId source,
        const T oldVals[],
        const PartitionId dest,
        const IndexType size ) const;

    template<typename T>
    void swapImpl( T val[], const IndexType n, PartitionId partner ) const;

    template<typename T>
    void maxlocImpl( T& val, int& location, PartitionId root ) const;

    template<typename T>
    void exchangeByPlanImpl(
        T* const recvData,
        const CommunicationPlan& recvPlan,
        const T* const sendData,
        const CommunicationPlan& sendPlan ) const;

    template<typename T>
    std::auto_ptr<SyncToken> exchangeByPlanAsyncImpl(
        T* const recvData,
        const CommunicationPlan& recvPlan,
        const T* const sendData,
        const CommunicationPlan& sendPlan ) const;

    const PGASInterface* mPGASInterface;
    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static const int defaultTag;

    virtual ContextPtr getCommunicationContext() const;
};
}

#endif // LAMA_PGASCOMMUNICATOR_H_
