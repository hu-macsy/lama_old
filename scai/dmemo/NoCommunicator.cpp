/**
 * @file NoCommunicator.cpp
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
 * @brief NoCommunicator.cpp
 * @author Thomas Brandes
 * @date 15.03.2011
 */

// hpp
#include <scai/dmemo/NoCommunicator.hpp>
#include <scai/dmemo/NoCollectiveFile.hpp>

// local library
#include <scai/dmemo/CommunicationPlan.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/safer_memcpy.hpp>

#include <cstring>
#include <memory>
#include <unistd.h>

using namespace std;

using scai::common::safer_memcpy;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( NoCommunicator::logger, "Communicator.NoCommunicator" )

NoCommunicator::NoCommunicator() : Communicator( CommunicatorType::NO )
{
    SCAI_LOG_DEBUG( logger, "NoCommunicator()" )
}

NoCommunicator::~NoCommunicator()
{
    SCAI_LOG_DEBUG( logger, "~NoCommunicator()" )
}

hmemo::ContextPtr NoCommunicator::getCommunicationContext( const hmemo::_HArray& ) const
{
    return hmemo::Context::getHostPtr();
}

std::unique_ptr<CollectiveFile> NoCommunicator::collectiveFile() const
{
    return std::unique_ptr<CollectiveFile>( new NoCollectiveFile( shared_from_this() ) );
}

bool NoCommunicator::isEqual( const Communicator& other ) const
{
    return typeid( *this ) == typeid( other );
}

Communicator::ThreadSafetyLevel NoCommunicator::getThreadSafetyLevel() const
{
    return Communicator::Multiple;
}

void NoCommunicator::all2allImpl( void* recvBuffer, const void* sendBuffer, const common::ScalarType stype ) const
{
    // exchange one single element on a single processor

    safer_memcpy( recvBuffer, sendBuffer, typeSize( stype ) );
}

void NoCommunicator::reduceImpl( 
    void* outData, 
    const void* inData, 
    const IndexType n, 
    const common::ScalarType stype, 
    const common::BinaryOp ) const
{
    if ( outData == inData )
    {
        return;   // IN_PLACE
    }

    safer_memcpy( outData, inData, typeSize( stype ) * n );
}

void NoCommunicator::scanImpl( void* outData, const void* inData, const IndexType n, const common::ScalarType stype ) const
{
    if ( outData == inData )
    {
        return;   // IN_PLACE
    }

    safer_memcpy( outData, inData, typeSize( stype ) * n );
}

/* ---------------------------------------------------------------------------------- */
/*      exchangeByPlan                                                                */
/* ---------------------------------------------------------------------------------- */

void NoCommunicator::exchangeByPlanImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    const common::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( recvPlan.size(), sendPlan.size(), "size mismatch" )

    if ( 0 == recvPlan.size() && 0 == sendPlan.size() )
    {
        return;
    }

    // send / recv plan have maximal one value
    SCAI_ASSERT_EQ_ERROR( 1, recvPlan.size(), "maximal one value in recvPlan" )

    IndexType quantity = recvPlan[0].quantity;

    // recv and send plan must have same quantity
    SCAI_ASSERT_EQ_ERROR( quantity, sendPlan[0].quantity, "quantity mismatch" )

    // self copy of send data to recv data

    safer_memcpy( recvData, sendData, quantity * common::typeSize( stype ) );
}

tasking::SyncToken* NoCommunicator::exchangeByPlanAsyncImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    const common::ScalarType stype ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, stype );
    return new tasking::NoSyncToken();
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

IndexType NoCommunicator::shiftImpl(
    void*,
    const IndexType,
    const PartitionId,
    const void*,
    const IndexType,
    const PartitionId,
    const common::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "shiftImpl should never be called for NoCommunicator" )
}

tasking::SyncToken* NoCommunicator::shiftAsyncImpl(
    void*,
    const PartitionId,
    const void*,
    const PartitionId,
    const IndexType,
    const common::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "shiftAsyncImpl should never be called for NoCommunicator" )
}

/* ---------------------------------------------------------------------------------- */
/*              minloc/maxloc                                                         */
/* ---------------------------------------------------------------------------------- */

void NoCommunicator::maxlocImpl( void*, IndexType*, PartitionId root, common::ScalarType ) const
{
    // nothing to do
    SCAI_ASSERT_EQ_ERROR( root, 0, "illegal root partition" )
}

void NoCommunicator::minlocImpl( void*, IndexType*, PartitionId root, common::ScalarType ) const
{
    // nothing to do
    SCAI_ASSERT_EQ_ERROR( root, 0, "illegal root partition" )
}

bool NoCommunicator::supportsLocReduction( common::ScalarType, common::ScalarType ) const
{
    return true;
}

/* ---------------------------------------------------------------------------------- */

void NoCommunicator::bcastImpl( void*, const IndexType, const PartitionId root, common::ScalarType ) const
{
    // Nothing to do as root is the only one processor
    SCAI_ASSERT_EQ_ERROR( root, 0, "root can only be processor 0" )
}

/* ---------------------------------------------------------------------------------- */

void NoCommunicator::all2allvImpl( 
    void* recvBuffer[], 
    const IndexType recvSizes[], 
    const void* sendBuffer[], 
    const IndexType sendSizes[], 
    const common::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( sendSizes[0], recvSizes[0], "serious mismatch" )

    SCAI_LOG_INFO( logger, "all2allv<" << stype << ">, size = " << sendSizes[0] )
    safer_memcpy( recvBuffer[0], sendBuffer[0], typeSize( stype ) * sendSizes[0] );
}

void NoCommunicator::scatterImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const common::ScalarType stype ) const
{
    SCAI_LOG_INFO( logger, "scatter<" << stype << ">, size = " << n << ", root = " << root )

    SCAI_ASSERT_EQ_ERROR( root, 0, "" )

    safer_memcpy( myVals, allVals, common::typeSize( stype ) * n );
}

void NoCommunicator::scatterVImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const IndexType sizes[],
    const common::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )
    SCAI_ASSERT_EQ_ERROR( sizes[0], n , "size mismatch" )

    safer_memcpy( myVals, allVals, common::typeSize( stype ) * n );
}

void NoCommunicator::gatherImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    common::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )

    safer_memcpy( allVals, myVals, common::typeSize( stype ) * n );
}

void NoCommunicator::gatherVImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const IndexType sizes[],
    common::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )
    SCAI_ASSERT_EQ_ERROR( sizes[0], n, "" )

    safer_memcpy( allVals, myVals, common::typeSize( stype ) * n );
}

void NoCommunicator::swapImpl( void*, const IndexType, const PartitionId partner, const common::ScalarType ) const
{
    SCAI_ASSERT_EQ_ERROR( partner, 0, "" )
}

/* --------------------------------------------------------------- */

void NoCommunicator::getProcessorName( char* name ) const
{
    size_t len = maxProcessorName();

    memset( name, 0, len * sizeof( char ) );

    gethostname( name, len );
}

size_t NoCommunicator::maxProcessorName() const
{
    return 256;
}

/* --------------------------------------------------------------- */

void NoCommunicator::synchronize() const
{
    // no synchronization needed
}

/* --------------------------------------------------------------- */

NoCommunicator* NoCommunicator::splitIt( PartitionId, PartitionId ) const
{
    return new NoCommunicator();
}

/* --------------------------------------------------------------- */

void NoCommunicator::writeAt( std::ostream& stream ) const
{
    stream << "NoComm";
}

/* --------------------------------------------------------------- */

static std::weak_ptr<class NoCommunicator> theNoCommunicatorInstance;

CommunicatorPtr NoCommunicator::create()
{
    std::shared_ptr<NoCommunicator> communicator;

    // use the last communicatorInstance if it is still valid

    if ( theNoCommunicatorInstance.expired() )
    {
        // create a new instance of NoCommunicator and keep it for further uses
        communicator = std::shared_ptr<NoCommunicator>( new NoCommunicator() );
        theNoCommunicatorInstance = communicator;
    }
    else
    {
        // the last communicator instance is still valid, so we return new shared pointer to it
        communicator = theNoCommunicatorInstance.lock();
    }

    return communicator;
}

/* --------------------------------------------------------------- */

CommunicatorType NoCommunicator::createValue()
{
    return CommunicatorType::NO;
}

} /* end namespace dmemo */

} /* end namespace scai */
