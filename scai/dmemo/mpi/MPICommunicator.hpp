/**
 * @file MPICommunicator.hpp
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
 * @brief MPICommunicator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 * @since 1.0.0
 */

#pragma once

#include <mpi.h> //Intel MPI need mpi.h to be included before stdio.h so this header comes first

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/dmemo/CRTPCommunicator.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Thread.hpp>
#include <scai/common/macros/unused.hpp>

// std
#include <vector>

namespace scai
{

namespace lama
{

/** Communicator class that implements communication and data exchange via MPI.
 *
 *  MPI_Init is called in the constructor, MPI_Finalize is called in the destructor.
 */

class COMMON_DLL_IMPORTEXPORT MPICommunicator: 

    public CRTPCommunicator<MPICommunicator>,
    public Communicator::Register<MPICommunicator>           // register at factory
{

// Only MPICommunicatorManager is allowed to create MPI communicator

    friend class MPICommunicatorManager;
    friend class CRTPCommunicator<MPICommunicator> ;

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
    virtual void all2all( IndexType recvValues[], const IndexType sendValues[] ) const;

    template<typename ValueType>
    MPI_Request startrecv( ValueType* buffer, int count, int source ) const;

    template<typename ValueType>
    MPI_Request startsend( const ValueType buffer[], int count, int target ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    MPICommunicator( int& argc, char** & argv );

    template<typename ValueType>
    inline static MPI_Datatype getMPIType();

    template<typename T1,typename T2>
    inline static MPI_Datatype getMPI2Type();

    template<typename ValueType>
	inline static MPI_Op getMPISum();

    template<typename ValueType>
    void bcastImpl( ValueType val[], const IndexType n, const PartitionId root ) const;

    template<typename ValueType>
    void all2allvImpl( ValueType* recvBuffer[], IndexType recvCount[], ValueType* sendBuffer[], IndexType sendCount[] ) const;

    template<typename ValueType>
    ValueType sumImpl( ValueType myVal ) const;

    template<typename ValueType>
    ValueType maxImpl( ValueType myVal ) const;

    template<typename ValueType>
    ValueType minImpl( ValueType myVal ) const;

    template<typename ValueType>
    inline void send( const ValueType* buffer, int count, int target ) const;

    template<typename ValueType>
    inline int getCount( MPI_Status& status ) const;

    template<typename ValueType>
    void scatterImpl( ValueType myvals[], const IndexType n, const PartitionId root, const ValueType allvals[] ) const;

    template<typename ValueType>
    void scatterVImpl(
        ValueType myvals[],
        const IndexType n,
        const PartitionId root,
        const ValueType allvals[],
        const IndexType sizes[] ) const;

    template<typename ValueType>
    void gatherImpl( ValueType allvals[], const IndexType n, const PartitionId root, const ValueType myvals[] ) const;

    template<typename ValueType>
    void gatherVImpl(
        ValueType allvals[],
        const IndexType n,
        const PartitionId root,
        const ValueType myvals[],
        const IndexType sizes[] ) const;

    template<typename ValueType>
    IndexType shiftImpl(
        ValueType newvals[],
        const IndexType newSize,
        const PartitionId source,
        const ValueType oldVals[],
        const IndexType oldSize,
        const PartitionId dest ) const;

    template<typename ValueType>
    tasking::SyncToken* shiftAsyncImpl(
        ValueType newvals[],
        const PartitionId source,
        const ValueType oldVals[],
        const PartitionId dest,
        const IndexType size ) const;

    template<typename ValueType>
    void swapImpl( ValueType val[], const IndexType n, PartitionId partner ) const;

    template<typename ValueType>
    void maxlocImpl( ValueType& val, IndexType& location, PartitionId root ) const;

    template<typename ValueType>
    void exchangeByPlanImpl(
        ValueType recvData[],
        const CommunicationPlan& recvPlan,
        const ValueType sendData[],
        const CommunicationPlan& sendPlan ) const;

    template<typename ValueType>
    tasking::SyncToken* exchangeByPlanAsyncImpl(
        ValueType recvData[],
        const CommunicationPlan& recvPlan,
        const ValueType sendData[],
        const CommunicationPlan& sendPlan ) const;

    void initialize( int& argc, char** & argv );

    const common::Thread::Id mMainThread;  // id of thread that calls constructor

    inline MPI_Comm selectMPIComm() const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    const int defaultTag;

protected:

    MPICommunicator( int& argc, char** & argv, const communicator::CommunicatorKind& type );

    MPICommunicator();

    void setNodeData();

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    int mRank; // rank of this processor
    int mSize;// size of communicator

    bool mExternInitialization;

    MPI_Comm mCommWorld;
    MPI_Comm mComm;
    MPI_Comm mCommTask;

    static MPI_Op mSumComplexLongDouble;

    static void sum_complex_long_double(void *in, void *out, int *count,
                                     MPI_Datatype *dtype);

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

    bool isCUDAAware;// if true data on CUDA context can be communicated

public:

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory 

    static communicator::CommunicatorKind createValue();

private:

    // Guard class whose destructor takes care of MPI finalize

    class MPIGuard
    {
        MPIGuard();
        ~MPIGuard();
    };

    static MPIGuard guard;   // define one guard variable
};

/* ---------------------------------------------------------------------------------- */
/*              getMPIType                                                            */
/* ---------------------------------------------------------------------------------- */

template<>
inline MPI_Datatype MPICommunicator::getMPIType<float>()
{
    return MPI_FLOAT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<char>()
{
    return MPI_CHAR;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<double>()
{
    return MPI_DOUBLE;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<int>()
{
    return MPI_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<unsigned int>()
{
    return MPI_UNSIGNED;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<long double>()
{
    return MPI_LONG_DOUBLE;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
inline MPI_Datatype MPICommunicator::getMPIType<ComplexFloat>()
{
    return MPI_COMPLEX;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<ComplexDouble>()
{
    // Be careful: some MPI implementations do not provide MPI_DOUBLE_COMPLEX
    // May be helpful: MPI::DOUBLE_COMPLEX
    return MPI_DOUBLE_COMPLEX;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<ComplexLongDouble>()
{
    // Be careful: some MPI implementations do not provide MPI_DOUBLE_COMPLEX
    // May be helpful: MPI::DOUBLE_COMPLEX
    return MPI::LONG_DOUBLE_COMPLEX;
}

#endif

template<>
inline MPI_Datatype MPICommunicator::getMPIType<unsigned long>()
{
    return MPI_UNSIGNED_LONG;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPI2Type                                                           */
/* ---------------------------------------------------------------------------------- */

template<typename T1,typename T2>
inline MPI_Datatype MPICommunicator::getMPI2Type()
{
    COMMON_THROWEXCEPTION( "unsupported type for MPI communication" )
    return MPI_2INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<float,int>()
{
    return MPI_FLOAT_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<double,int>()
{
    return MPI_DOUBLE_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<int,int>()
{
    return MPI_2INT;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPISum                                                           */
/* ---------------------------------------------------------------------------------- */


template<typename ValueType>
inline MPI_Op MPICommunicator::getMPISum()
{
    return MPI_SUM;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
inline MPI_Op MPICommunicator::getMPISum<ComplexLongDouble>()
{
    return mSumComplexLongDouble;
}

#endif

} /* end namespace lama */

} /* end namespace scai */
