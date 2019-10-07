/**
 * @file NoCommunicator.hpp
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
 * @brief Implementation of a communicator class for non-distributed (replicated) objects.
 * @author Thomas Brandes
 * @date 15.03.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Communicator.hpp>

namespace scai
{

namespace dmemo
{

/** 
 *  The class NoCommunicator stands for a communicator that contains just this single processor.
 *
 *  It is used a fallback communicator if no MPI is supported, but is also used in cases
 *  where trivial communication is required.
 */

class COMMON_DLL_IMPORTEXPORT NoCommunicator:

    public Communicator,
    public Communicator::Register<NoCommunicator>           // register at factory
{
public:

    NoCommunicator();

    virtual ~NoCommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    /** Implementation of pure method Communciator::all2allImpl */

    virtual void all2allImpl( void* recvBuffer, const void* sendBuffer, const common::ScalarType stype ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

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

    /** Implementation of pure method Communicator::bcastImpl */

    void bcastImpl( void* val, const IndexType n, const PartitionId root, const common::ScalarType stype ) const;

    /** Implementation of pure method Communciator::all2allvImpl */

    void all2allvImpl( void* recvBuffer[], const IndexType recvCount[],
                       const void* sendBuffer[], const IndexType sendCount[],
                       const common::ScalarType stype ) const;

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

    void maxlocImpl( void* val, IndexType* location, PartitionId root, common::ScalarType stype ) const;

    void minlocImpl( void* val, IndexType* location, PartitionId root, common::ScalarType stype ) const;

    /** Implementation of Communicator::supportsLocReduction */

    virtual bool supportsLocReduction( common::ScalarType vType, common::ScalarType iType ) const;

    void swapImpl( void* val, const IndexType n, const PartitionId partner, common::ScalarType stype ) const;

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

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    virtual std::unique_ptr<class CollectiveFile> collectiveFile() const;

    /** Implementation of Communicator::reduceImpl */

    virtual void reduceImpl( 
        void* outValues, 
        const void* inValues, 
        const IndexType n, 
        const common::ScalarType stype,
        const common::BinaryOp op ) const;

    /** Implementation of Communicator::scanImpl */

    virtual void scanImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const;

    /** Implementation of Communicator::getProcessorName */

    virtual void getProcessorName( char* name ) const;

    /** Implementation of Communicator::maxProcessorName */

    virtual size_t maxProcessorName() const;

public:

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory

    static CommunicatorType createValue();

protected:

    /** Implementation of pure method NoCommunicator::splitIt */

    virtual NoCommunicator* splitIt( PartitionId color, PartitionId key ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace dmemo */

} /* end namespace scai */
