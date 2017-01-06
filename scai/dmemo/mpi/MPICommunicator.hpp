/**
 * @file MPICommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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

    public Communicator,
    public Communicator::Register<MPICommunicator>           // register at factory
{
public:

    virtual ~MPICommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

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

    MPI_Request startrecv( void* buffer, int count, int source, common::scalar::ScalarType stype ) const;

    MPI_Request startsend( const void* buffer, int count, int target, common::scalar::ScalarType stype ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    MPICommunicator( int& argc, char** & argv );

    /** Implementation of Communicator::sumImpl */

    virtual void sumImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::minImpl */

    virtual void minImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::maxImpl */

    virtual void maxImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Translate SCAI project enum type ScalarType to MPI enum MPI_Datatype */

    inline static MPI_Datatype getMPIType( common::scalar::ScalarType stype );

    /** Translate SCAI project enum type ScalarType ( two types ) to MPI enum MPI_Datatype */

    inline static MPI_Datatype getMPI2Type( common::scalar::ScalarType stype1, common::scalar::ScalarType stype2 );

    inline static MPI_Op getMPISum( common::scalar::ScalarType stype );

    inline static MPI_Op getMPIMin( common::scalar::ScalarType stype );

    inline static MPI_Op getMPIMax( common::scalar::ScalarType stype );

    /** MPI implementation for pure method Communicator::bcastImpl */

    void bcastImpl( void* val, const IndexType n, const PartitionId root, common::scalar::ScalarType stype ) const;

    /** MPI Implementation for pure method Communciator::all2allvImpl */

    void all2allvImpl( void* recvBuffer[], IndexType recvCount[],
                       void* sendBuffer[], IndexType sendCount[],
                       common::scalar::ScalarType stype ) const;

    inline void send( const void* buffer, int count, int target, common::scalar::ScalarType stype ) const;

    inline int getCount( MPI_Status& status, common::scalar::ScalarType stype ) const;

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

    void maxlocImpl( void* val, IndexType* location, PartitionId root, common::scalar::ScalarType stype ) const;

    void minlocImpl( void* val, IndexType* location, PartitionId root, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::supportsLocReduction */

    virtual bool supportsLocReduction( common::scalar::ScalarType vType, common::scalar::ScalarType iType ) const;

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

    /** Implementation of Communicator::getProcessorName */

    virtual void getProcessorName( char* name ) const;

    /** Implementation of Communicator::maxProcessorName */

    virtual size_t maxProcessorName() const;

    void initialize( int& argc, char** & argv );

    const common::Thread::Id mMainThread;  // id of thread that calls constructor

    inline MPI_Comm selectMPIComm() const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    const int defaultTag;

protected:

    MPICommunicator( int& argc, char** & argv, const CommunicatorKind& type );

    MPICommunicator();

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    bool mExternInitialization;

    MPI_Comm mCommWorld;
    MPI_Comm mComm;
    MPI_Comm mCommTask;

#ifdef SCAI_COMPLEX_SUPPORTED
    static MPI_Op mSumComplexLongDouble;

    static void sum_complex_long_double( void* in, void* out, int* count,
                                         MPI_Datatype* dtype );

    static MPI_Datatype mComplexLongDoubleType;

    static MPI_Op mMaxComplexFloat;
    static MPI_Op mMaxComplexDouble;
    static MPI_Op mMaxComplexLongDouble;

    template<typename ValueType>
    static void max_operator( void* in, void* out, int* count, MPI_Datatype* dtype );

    static MPI_Op mMinComplexFloat;
    static MPI_Op mMinComplexDouble;
    static MPI_Op mMinComplexLongDouble;

    template<typename ValueType>
    static void min_operator( void* in, void* out, int* count, MPI_Datatype* dtype );
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

inline MPI_Datatype MPICommunicator::getMPIType( common::scalar::ScalarType stype )
{
    switch ( stype )
    {
        case common::scalar::INT                 :
            return MPI_INT;
        case common::scalar::LONG                :
            return MPI_LONG;
        case common::scalar::FLOAT               :
            return MPI_FLOAT;
        case common::scalar::DOUBLE              :
            return MPI_DOUBLE;
        case common::scalar::LONG_DOUBLE         :
            return MPI_LONG_DOUBLE;
        case common::scalar::COMPLEX             :
            return MPI_COMPLEX;
        case common::scalar::DOUBLE_COMPLEX      :
            return MPI_DOUBLE_COMPLEX;
        case common::scalar::LONG_DOUBLE_COMPLEX :
            return mComplexLongDoubleType;
        case common::scalar::CHAR                :
            return MPI_CHAR;
        case common::scalar::UNSIGNED_INT        :
            return MPI_UNSIGNED;
        case common::scalar::UNSIGNED_LONG       :
            return MPI_UNSIGNED_LONG;

        default:
            COMMON_THROWEXCEPTION( "No MPI Type specified for " << stype )
            return MPI_INT;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPI2Type                                                           */
/* ---------------------------------------------------------------------------------- */

inline MPI_Datatype MPICommunicator::getMPI2Type( common::scalar::ScalarType stype1,
        common::scalar::ScalarType stype2 )
{
    if ( stype2 != common::scalar::INT )
    {
        COMMON_THROWEXCEPTION( "getMPI2Type, 2nd type must be INT, is " << stype2 )
    }

    switch ( stype1 )
    {
        case common::scalar::INT                 :
            return MPI_2INT;
        case common::scalar::FLOAT               :
            return MPI_FLOAT_INT;
        case common::scalar::DOUBLE              :
            return MPI_DOUBLE_INT;

        default:
            COMMON_THROWEXCEPTION( "No MPI2 Type for " << stype1 )
            return MPI_2INT;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPISum                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPISum( common::scalar::ScalarType stype )
{

#ifdef SCAI_COMPLEX_SUPPORTED

    if ( stype == common::scalar::LONG_DOUBLE_COMPLEX )
    {
        return mSumComplexLongDouble;
    }

#endif

    return MPI_SUM;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPIMax                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPIMax( common::scalar::ScalarType stype )
{
    switch ( stype )
    {
#ifdef SCAI_COMPLEX_SUPPORTED
        case common::scalar::COMPLEX             :
            return mMaxComplexFloat;
        case common::scalar::DOUBLE_COMPLEX      :
            return mMaxComplexDouble;
        case common::scalar::LONG_DOUBLE_COMPLEX :
            return mMaxComplexLongDouble;
#endif
        default:
            return MPI_MAX;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPIMin                                                             */
/* ---------------------------------------------------------------------------------- */

inline MPI_Op MPICommunicator::getMPIMin( common::scalar::ScalarType stype )
{
    switch ( stype )
    {
#ifdef SCAI_COMPLEX_SUPPORTED
        case common::scalar::COMPLEX             :
            return mMinComplexFloat;
        case common::scalar::DOUBLE_COMPLEX      :
            return mMinComplexDouble;
        case common::scalar::LONG_DOUBLE_COMPLEX :
            return mMinComplexLongDouble;
#endif
        default:
            return MPI_MIN;
    }
}

} /* end namespace dmemo */

} /* end namespace scai */
