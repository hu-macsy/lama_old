/**
 * @file MPICommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
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
#include <scai/dmemo/Communicator.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/unused.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

enum class MPICommKind
{
    EXTERNAL,   //!<  communicator was already initialized 
    INTERNAL,   //!<  comm world created
    CREATED     //!<  comm created
};

/** Communicator class that implements communication and data exchange via MPI.
 *
 *  MPI_Init is called in the constructor, MPI_Finalize is called in the destructor.
 */

class COMMON_DLL_IMPORTEXPORT MPICommunicator:

    public Communicator,
    public Communicator::Register<MPICommunicator>           // register at factory
{
public:

    virtual ~MPICommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    MPI_Comm getMPIComm() const;

    /** MPI Implementation for pure method Communciator::all2allImpl */

    void all2allImpl( void* recvBuffer, 
                      const void* sendBuffer,
                      const common::ScalarType stype ) const;

    MPI_Request startrecv( void* buffer, int count, int source, const common::ScalarType stype ) const;

    MPI_Request startsend( const void* buffer, int count, int target, const common::ScalarType stype ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    MPICommunicator( int& argc, char** & argv );

    /** Implementation of Communicator::sumImpl */

    virtual void sumImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const;

    /** Implementation of Communicator::minImpl */

    virtual void minImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const;

    /** Implementation of Communicator::maxImpl */

    virtual void maxImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const;

    /** Implementation of Communicator::scanImpl */

    virtual void scanImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const;

    /** Translate SCAI project enum type ScalarType to MPI enum MPI_Datatype */

    inline static MPI_Datatype getMPIType( const common::ScalarType stype );

    /** Translate SCAI project enum type ScalarType ( two types ) to MPI enum MPI_Datatype */

    inline static MPI_Datatype getMPI2Type( const common::ScalarType stype1, const common::ScalarType stype2 );

    inline static MPI_Op getMPISum( const common::ScalarType stype );

    inline static MPI_Op getMPIMin( const common::ScalarType stype );

    inline static MPI_Op getMPIMax( const common::ScalarType stype );

    /** MPI implementation for pure method Communicator::bcastImpl */

    void bcastImpl( void* val, const IndexType n, const PartitionId root, const common::ScalarType stype ) const;

    /** MPI Implementation for pure method Communciator::all2allvImpl */

    void all2allvImpl( void* recvBuffer[], const IndexType recvCount[],
                       const void* sendBuffer[], const IndexType sendCount[],
                       const common::ScalarType stype ) const;

    inline void send( const void* buffer, int count, int target, const common::ScalarType stype ) const;

    inline int getCount( MPI_Status& status, const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::scatterImpl */

    void scatterImpl( void* myVals, const IndexType n, const PartitionId root, const void* allVals, const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::scatterVImpl */

    void scatterVImpl( void* myVals, const IndexType n, const PartitionId root,
                       const void* allVals, const IndexType sizes[], const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherImpl */

    void gatherImpl( void* allVals, const IndexType n, const PartitionId root, const void* myVals, const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherVImpl */

    void gatherVImpl(
        void* allvals,
        const IndexType n,
        const PartitionId root,
        const void* myvals,
        const IndexType sizes[],
        const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::shiftImpl */

    IndexType shiftImpl(
        void* newVals,
        const IndexType newSize,
        const PartitionId source,
        const void* oldVals,
        const IndexType oldSize,
        const PartitionId dest,
        const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::shiftAsyncImpl */

    tasking::SyncToken* shiftAsyncImpl(
        void* newVals,
        const PartitionId source,
        const void* oldVals,
        const PartitionId dest,
        const IndexType size,
        const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::swapImpl */

    void swapImpl( void* val, const IndexType n, PartitionId partner, const common::ScalarType stype ) const;

    void maxlocImpl( void* val, IndexType* location, PartitionId root, const common::ScalarType stype ) const;

    void minlocImpl( void* val, IndexType* location, PartitionId root, const common::ScalarType stype ) const;

    /** Implementation of Communicator::supportsLocReduction */

    virtual bool supportsLocReduction( const common::ScalarType vType, const common::ScalarType iType ) const;

    /** Implementation of pure method Communicator::exchangeByPlanImpl */

    void exchangeByPlanImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::ScalarType stype ) const;

    /** Implementation of pure method Communicator::exchangeByPlanAsyncImpl */

    tasking::SyncToken* exchangeByPlanAsyncImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::ScalarType stype ) const;

    /** Implementation of Communicator::getProcessorName */

    virtual void getProcessorName( char* name ) const;

    /** Implementation of Communicator::maxProcessorName */

    virtual size_t maxProcessorName() const;

    void initialize( int& argc, char** & argv );

    const std::thread::id mMainThread;  // id of thread that calls constructor

    inline MPI_Comm selectMPIComm() const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    const int defaultTag;

protected:

    /** Implementation of pure method Communicator::splitIt */

    virtual MPICommunicator* splitIt( PartitionId color, PartitionId key ) const;

    MPICommunicator( int& argc, char** & argv, const CommunicatorKind& type );

    MPICommunicator();

    MPICommunicator( const MPICommunicator& comm, const PartitionId color, const PartitionId key );

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    MPICommKind mKind;    // kind of communicator needed for destructor

    MPI_Comm mComm;

#ifdef SCAI_COMPLEX_SUPPORTED
    static MPI_Op mSumComplexLongDouble;

    static void sum_complex_long_double( void* in, void* out, int* count,
                                         MPI_Datatype* dtype );

    static MPI_Datatype mComplexLongDoubleType;

#endif

    Communicator::ThreadSafetyLevel mThreadSafetyLevel;

    bool isCUDAAware;   // if true data on CUDA context can be communicated

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

inline MPI_Datatype MPICommunicator::getMPIType( const common::ScalarType stype )
{
    switch ( stype )
    {
        case common::ScalarType::INT                 :
            return MPI_INT;
        case common::ScalarType::LONG                :
            return MPI_LONG;
        case common::ScalarType::FLOAT               :
            return MPI_FLOAT;
        case common::ScalarType::DOUBLE              :
            return MPI_DOUBLE;
        case common::ScalarType::LONG_DOUBLE         :
            return MPI_LONG_DOUBLE;
#ifdef SCAI_COMPLEX_SUPPORTED

#ifndef MPI_C_COMPLEX
   #define MPI_C_COMPLEX MPI_COMPLEX
#endif

#ifndef MPI_C_DOUBLE_COMPLEX
   #define MPI_C_DOUBLE_COMPLEX MPI_DOUBLE_COMPLEX
#endif

        case common::ScalarType::COMPLEX             :
            return MPI_C_COMPLEX;
        case common::ScalarType::DOUBLE_COMPLEX      :
            return MPI_C_DOUBLE_COMPLEX;
        case common::ScalarType::LONG_DOUBLE_COMPLEX :
            return mComplexLongDoubleType;
#endif
        case common::ScalarType::CHAR                :
            return MPI_CHAR;
        case common::ScalarType::UNSIGNED_INT        :
            return MPI_UNSIGNED;
        case common::ScalarType::UNSIGNED_LONG       :
            return MPI_UNSIGNED_LONG;

        default:
            COMMON_THROWEXCEPTION( "No MPI Type specified for " << stype )
            return MPI_INT;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPI2Type                                                           */
/* ---------------------------------------------------------------------------------- */

inline MPI_Datatype MPICommunicator::getMPI2Type( 
    const common::ScalarType stype1,
    const common::ScalarType stype2 )
{
    if ( stype2 != common::ScalarType::INT )
    {
        COMMON_THROWEXCEPTION( "getMPI2Type, 2nd type must be INT, is " << stype2 )
    }

    switch ( stype1 )
    {
        case common::ScalarType::INT                 :
            return MPI_2INT;
        case common::ScalarType::FLOAT               :
            return MPI_FLOAT_INT;
        case common::ScalarType::DOUBLE              :
            return MPI_DOUBLE_INT;

        default:
            COMMON_THROWEXCEPTION( "No MPI2 Type for " << stype1 )
            return MPI_2INT;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPISum                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPISum( const common::ScalarType stype )
{
    if ( stype == common::ScalarType::LONG_DOUBLE_COMPLEX )
    {

#ifdef SCAI_COMPLEX_SUPPORTED
        return mSumComplexLongDouble;
#else
        COMMON_THROWEXCEPTION( "MPI sum unsupported for " << stype )
        return MPI_SUM;
#endif
    }

    return MPI_SUM;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPIMax                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPIMax( const common::ScalarType stype )
{
    if ( common::isComplex( stype ) )
    {
        COMMON_THROWEXCEPTION( "max for complex types unsupported" )
    }

    return MPI_MAX;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPIMin                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPIMin( const common::ScalarType stype )
{
    if ( common::isComplex( stype ) )
    {
        COMMON_THROWEXCEPTION( "max for complex types unsupported" )
    }

    return MPI_MIN;
}

} /* end namespace dmemo */

} /* end namespace scai */
