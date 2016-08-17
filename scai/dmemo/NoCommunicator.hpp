/**
 * @file NoCommunicator.hpp
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
 * @brief Implementation of a communicator class for non-distributed (replicated) objects.
 * @author Thomas Brandes
 * @date 15.03.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/CRTPCommunicator.hpp>

namespace scai
{

namespace dmemo
{

/** The class NoCommunicator stands for objects that are replicated on each
 *  partition or processor.
 */

class COMMON_DLL_IMPORTEXPORT NoCommunicator:

    public CRTPCommunicator<NoCommunicator>,
    public Communicator::Register<NoCommunicator>           // register at factory
{

    friend class CRTPCommunicator<NoCommunicator> ;

public:

    NoCommunicator();

    virtual ~NoCommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    virtual PartitionId getSize() const;

    virtual PartitionId getRank() const;

    virtual PartitionId getNodeSize() const;

    virtual PartitionId getNodeRank() const;

    virtual void all2all( IndexType recvValues[], const IndexType sendValues[] ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    /** Implementation of pure method Communicator::shiftData */

    IndexType shiftData(
        void* newVals,
        const IndexType newSize,
        const PartitionId source,
        const void* oldVals,
        const IndexType oldSize,
        const PartitionId dest,
        common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::shiftAsyncData */

    tasking::SyncToken* shiftAsyncData(
        void* newVals,
        const PartitionId source,
        const void* oldVals,
        const PartitionId dest,
        const IndexType size,
        common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::bcastData */

    void bcastData( void* val, const IndexType n, const PartitionId root, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communciator::all2allvImpl */

    void all2allvImpl( void* recvBuffer[], IndexType recvCount[],
                       void* sendBuffer[], IndexType sendCount[],
                       common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::scatterData */

    void scatterData( void* myVals, const IndexType n, const PartitionId root, const void* allVals, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::scatterVData */

    void scatterVData( void* myVals, const IndexType n, const PartitionId root,
                       const void* allVals, const IndexType sizes[], common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherData */

    void gatherData( void* allVals, const IndexType n, const PartitionId root, const void* myVals, common::scalar::ScalarType stype ) const;

    /** Implementation of pure method Communicator::gatherVData */

    void gatherVData(
        void* allvals,
        const IndexType n,
        const PartitionId root,
        const void* myvals,
        const IndexType sizes[],
        const common::scalar::ScalarType stype ) const;

    template<typename T>
    T maxImpl( const T value ) const;

    template<typename T>
    T minImpl( const T value ) const;

    template<typename T>
    void maxlocImpl( T& val, IndexType& location, const PartitionId root ) const;

    void swapImpl( void* val, const IndexType n, const PartitionId partner, common::scalar::ScalarType stype ) const;

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

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

    /** Implementation of Communicator::sumData */

    virtual void sumData( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::minData */

    virtual void minData( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

    /** Implementation of Communicator::maxData */

    virtual void maxData( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const;

public:

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory

    static CommunicatorKind createValue();
};

} /* end namespace dmemo */

} /* end namespace scai */
