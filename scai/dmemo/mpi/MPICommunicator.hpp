/**
 * @file MPICommunicator.hpp
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
 * @endlicense
 *
 * @brief MPICommunicator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
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

namespace dmemo
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

    template<typename T1, typename T2>
    inline static MPI_Datatype getMPI2Type();

    template<typename ValueType>
    inline static MPI_Op getMPISum();

    template<typename ValueType>
    inline static MPI_Op getMPIMin();

    template<typename ValueType>
    inline static MPI_Op getMPIMax();

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

    MPICommunicator( int& argc, char** & argv, const CommunicatorKind& type );

    MPICommunicator();

    void setNodeData();

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    int mRank; // rank of this processor
    int mSize;// size of communicator

    bool mExternInitialization;

    MPI_Comm mCommWorld;
    MPI_Comm mComm;
    MPI_Comm mCommTask;

#ifdef SCAI_COMPLEX_SUPPORTED
    static MPI_Op mSumComplexLongDouble;

    static void sum_complex_long_double( void* in, void* out, int* count,
                                         MPI_Datatype* dtype );

    static MPI_Op mMaxComplexFloat;
    static MPI_Op mMaxComplexDouble;
    static MPI_Op mMaxComplexLongDouble;

    template<typename ValueType>
    static void max_operator( void* in, void *out, int *count, MPI_Datatype *dtype );

    static MPI_Op mMinComplexFloat;
    static MPI_Op mMinComplexDouble;
    static MPI_Op mMinComplexLongDouble;

    template<typename ValueType>
    static void min_operator( void* in, void *out, int *count, MPI_Datatype *dtype );
#endif

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

    bool isCUDAAware;// if true data on CUDA context can be communicated

public:

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory

    static CommunicatorKind createValue();

private:

    // Guard class whose destructor takes care of MPI finalize

    class MPIGuard
    {
    public:

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

template<typename T1, typename T2>
inline MPI_Datatype MPICommunicator::getMPI2Type()
{
    COMMON_THROWEXCEPTION( "unsupported type for MPI communication" )
    return MPI_2INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<float, int>()
{
    return MPI_FLOAT_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<double, int>()
{
    return MPI_DOUBLE_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<int, int>()
{
    return MPI_2INT;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPISum                                                             */
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

/* ---------------------------------------------------------------------------------- */
/*              getMPIMax                                                             */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
inline MPI_Op MPICommunicator::getMPIMax()
{
    return MPI_MAX;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
inline MPI_Op MPICommunicator::getMPIMax<ComplexFloat>()
{
    return mMaxComplexFloat;
}

template<>
inline MPI_Op MPICommunicator::getMPIMax<ComplexDouble>()
{
    return mMaxComplexDouble;
}

template<>
inline MPI_Op MPICommunicator::getMPIMax<ComplexLongDouble>()
{
    return mMaxComplexLongDouble;
}

#endif

/* ---------------------------------------------------------------------------------- */
/*              getMPIMin                                                             */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
inline MPI_Op MPICommunicator::getMPIMin()
{
    return MPI_MIN;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
inline MPI_Op MPICommunicator::getMPIMin<ComplexFloat>()
{
    return mMinComplexFloat;
}

template<>
inline MPI_Op MPICommunicator::getMPIMin<ComplexDouble>()
{
    return mMinComplexDouble;
}

template<>
inline MPI_Op MPICommunicator::getMPIMin<ComplexLongDouble>()
{
    return mMinComplexLongDouble;
}

#endif

} /* end namespace dmemo */

} /* end namespace scai */
