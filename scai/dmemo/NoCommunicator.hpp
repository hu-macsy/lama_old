/**
 * @file NoCommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

    // Implementation methods are all private, but CRTPCommunicator is a friend class

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
    void bcastImpl( T val[], const IndexType n, const PartitionId root ) const;

    template<typename ValueType>
    void all2allvImpl( ValueType* recvBuffer[], IndexType recvCount[], ValueType* sendBuffer[], IndexType sendCount[] ) const;

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
    T sumImpl( const T value ) const;

    template<typename T>
    T maxImpl( const T value ) const;

    template<typename T>
    T minImpl( const T value ) const;

    template<typename T>
    void maxlocImpl( T& val, IndexType& location, const PartitionId root ) const;

    template<typename T>
    void swapImpl( T val[], const IndexType n, const PartitionId partner ) const;

    // common implementation for self exchange, uses size of datatype

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

    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const;

public:

    // static methods, variables to register create routine in Communicator factory of base class.

    static CommunicatorPtr create();

    // key for factory

    static CommunicatorKind createValue();
};

} /* end namespace dmemo */

} /* end namespace scai */
