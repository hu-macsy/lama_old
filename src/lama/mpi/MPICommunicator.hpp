/**
 * @file MPICommunicator.hpp
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
 * @brief MPICommunicator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 * $Id$
 */

#ifndef LAMA_MPI_COMMUNICATOR_HPP_
#define LAMA_MPI_COMMUNICATOR_HPP_

#include <mpi.h> //Intel MPI need mpi.h to be included before stdio.h so this header comes first

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Communicator.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/thread/thread.hpp>

#include <vector>

namespace lama
{

/** Communicator class that implements communication and data exchange via MPI.
 *
 *  MPI_Init is called in the constructor, MPI_Finalize is called in the destructor.
 */

class LAMA_DLL_IMPORTEXPORT MPICommunicator: public Communicator
{

// Only MPICommunicatorManager is allowed to create MPI communicator

    friend class MPICommunicatorManager;

public:
    virtual ~MPICommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    virtual PartitionId getSize() const;

    virtual PartitionId getRank() const;

    MPI_Comm getMPIComm() const;

    /** All-to-all exchange of an integer value between all processors.
     *
     * @param[out] recvValues will contain one value from each processor
     * @param[in]  sendValues must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( int* recvValues, const int* sendValues ) const;

    /** Exchange data between all processors by communication plans.
     *
     *  @param[out] recvData  buffer for data received from other processors
     *  @param[in]  recvPlan  number of elements and offsets for receiving
     *  @param[in]  sendData  buffer for data to send to other processors
     *  @param[in]  sendPlan  contains number of elements and offsets for sending data
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

    template<typename T>
    MPI_Request startrecv( T* buffer, int count, int source ) const;

    template<typename T>
    MPI_Request startsend( const T* buffer, int count, int target ) const;

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
     *  @param[in,out] val      is a value on each processor, only out for root with maximal value
     *  @param[in,out] location is an additional int value, only out for root
     *  @param[in]     root     TODO[doxy] Complete Description.
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

    MPICommunicator( int& argc, char** & argv );

    template<typename T>
    inline static MPI_Datatype getMPIType();

    template<typename T1,typename T2>
    inline static MPI_Datatype getMPI2Type();

    template<typename T>
    T sumImpl( T myVal ) const;

    template<typename T>
    T maxval( T myVal ) const;

    template<typename T>
    T minval( T myVal ) const;

    template<typename T>
    inline void send( const T* buffer, int count, int target ) const;

    template<typename T>
    inline int getCount( MPI_Status& status ) const;

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

    MPI_Comm mCommWorld;
    MPI_Comm mComm;
    MPI_Comm mCommTask;

    const boost::thread mMainThread;

    inline MPI_Comm selectMPIComm() const;

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

    int mRank; // rank of this processor
    int mSize; // size of communicator

    bool mExternInitialization;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static const int defaultTag;
protected:
    virtual ContextPtr getCommunicationContext() const;

};

}

#endif // LAMA_MPI_COMMUNICATOR_HPP_
