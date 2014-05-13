/**
 * @file MPICommunicator.hpp
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
 * @brief MPICommunicator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 * @since 1.0.0
 */

#ifndef LAMA_MPI_COMMUNICATOR_HPP_
#define LAMA_MPI_COMMUNICATOR_HPP_

#include <mpi.h> //Intel MPI need mpi.h to be included before stdio.h so this header comes first

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/CRTPCommunicator.hpp>

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

class LAMA_DLL_IMPORTEXPORT MPICommunicator: public CRTPCommunicator<MPICommunicator>
{

// Only MPICommunicatorManager is allowed to create MPI communicator

    friend class MPICommunicatorManager;
    friend class CRTPCommunicator<MPICommunicator>;

public:
    virtual ~MPICommunicator();

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

    MPI_Comm getMPIComm() const;

    /** All-to-all exchange of an integer value between all processors.
     *
     * @param[out] recvValues will contain one value from each processor
     * @param[in]  sendValues must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( int recvValues[], const int sendValues[] ) const;

    template<typename T>
    MPI_Request startrecv( T* buffer, int count, int source ) const;

    template<typename T>
    MPI_Request startsend( const T buffer[], int count, int target ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    MPICommunicator( int& argc, char** & argv );

    template<typename T>
    inline static MPI_Datatype getMPIType();

    template<typename T1,typename T2>
    inline static MPI_Datatype getMPI2Type();

    template<typename T>
    void bcastImpl( T val[], const IndexType n, const PartitionId root ) const;

    template<typename T>
    T sumImpl( T myVal ) const;

    template<typename T>
    T maxImpl( T myVal ) const;

    template<typename T>
    T minImpl( T myVal ) const;

    template<typename T>
    inline void send( const T* buffer, int count, int target ) const;

    template<typename T>
    inline int getCount( MPI_Status& status ) const;

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
    SyncToken* shiftAsyncImpl(
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
        T recvData[],
        const CommunicationPlan& recvPlan,
        const T sendData[],
        const CommunicationPlan& sendPlan ) const;

    template<typename T>
    SyncToken* exchangeByPlanAsyncImpl(
        T recvData[],
        const CommunicationPlan& recvPlan,
        const T sendData[],
        const CommunicationPlan& sendPlan ) const;

    void initialize( int& argc, char** & argv );

    const boost::thread mMainThread;

    inline MPI_Comm selectMPIComm() const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static const int defaultTag;
protected:
    MPICommunicator( int& argc, char** & argv, const std::string& type );

    void setNodeData();

    virtual ContextPtr getCommunicationContext() const;

    int mRank; // rank of this processor
    int mSize; // size of communicator

    bool mExternInitialization;

    MPI_Comm mCommWorld;
    MPI_Comm mComm;
    MPI_Comm mCommTask;

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;
};

}

#endif // LAMA_MPI_COMMUNICATOR_HPP_
