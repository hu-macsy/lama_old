/**
 * @file NoCommunicator.cpp
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
 * @brief NoCommunicator.cpp
 * @author Thomas Brandes
 * @date 15.03.2011
 */

// hpp
#include <scai/dmemo/NoCommunicator.hpp>

// local library
#include <scai/dmemo/CommunicationPlan.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/weak_ptr.hpp>

#include <cstring>
#include <unistd.h>

using namespace std;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( NoCommunicator::logger, "Communicator.NoCommunicator" )

NoCommunicator::NoCommunicator() : Communicator( NO )
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

bool NoCommunicator::isEqual( const Communicator& other ) const
{
    return typeid( *this ) == typeid( other );
}

Communicator::ThreadSafetyLevel NoCommunicator::getThreadSafetyLevel() const
{
    return Communicator::Multiple;
}

void NoCommunicator::all2all( IndexType recvValues[], const IndexType sendValues[] ) const
{
    recvValues[0] = sendValues[0];
}

void NoCommunicator::sumImpl( void* outData, const void* inData, const IndexType n, common::scalar::ScalarType stype ) const
{
    if ( outData == inData )
    {
        return;   // IN_PLACE
    }

    memcpy( outData, inData, typeSize( stype ) * n );
}

void NoCommunicator::minImpl( void* outData, const void* inData, const IndexType n, common::scalar::ScalarType stype ) const
{
    if ( outData == inData )
    {
        return;   // IN_PLACE
    }

    memcpy( outData, inData, typeSize( stype ) * n );
}

void NoCommunicator::maxImpl( void* outData, const void* inData, const IndexType n, common::scalar::ScalarType stype ) const
{
    if ( outData == inData )
    {
        return;   // IN_PLACE
    }

    memcpy( outData, inData, typeSize( stype ) * n );
}

/* ---------------------------------------------------------------------------------- */
/*      exchangeByPlan                                                                */
/* ---------------------------------------------------------------------------------- */

void NoCommunicator::exchangeByPlanImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    const common::scalar::ScalarType stype ) const
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

    memcpy( recvData, sendData, quantity * common::typeSize( stype ) );
}

tasking::SyncToken* NoCommunicator::exchangeByPlanAsyncImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    const common::scalar::ScalarType stype ) const
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
    const common::scalar::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "shiftImpl should never be called for NoCommunicator" )
}

tasking::SyncToken* NoCommunicator::shiftAsyncImpl(
    void*,
    const PartitionId,
    const void*,
    const PartitionId,
    const IndexType,
    const common::scalar::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "shiftAsyncImpl should never be called for NoCommunicator" )
}

/* ---------------------------------------------------------------------------------- */
/*              minloc/maxloc                                                         */
/* ---------------------------------------------------------------------------------- */

void NoCommunicator::maxlocImpl( void*, IndexType*, PartitionId root, common::scalar::ScalarType ) const
{
    // nothing to do
    SCAI_ASSERT_EQ_ERROR( root, 0, "illegal root partition" )
}

void NoCommunicator::minlocImpl( void*, IndexType*, PartitionId root, common::scalar::ScalarType ) const
{
    // nothing to do
    SCAI_ASSERT_EQ_ERROR( root, 0, "illegal root partition" )
}

bool NoCommunicator::supportsLocReduction( common::scalar::ScalarType, common::scalar::ScalarType ) const
{
    return true;  
}

/* ---------------------------------------------------------------------------------- */

void NoCommunicator::bcastImpl( void*, const IndexType, const PartitionId root, common::scalar::ScalarType ) const
{
    // Nothing to do as root is the only one processor
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )
}

/* ---------------------------------------------------------------------------------- */

void NoCommunicator::all2allvImpl( void*[], IndexType[], void*[], IndexType[], common::scalar::ScalarType ) const
{
}

void NoCommunicator::scatterImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const common::scalar::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )

    memcpy( myVals, allVals, common::typeSize( stype ) * n );
}

void NoCommunicator::scatterVImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const IndexType sizes[],
    const common::scalar::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )
    SCAI_ASSERT_EQ_ERROR( sizes[0], n , "size mismatch" )

    memcpy( myVals, allVals, common::typeSize( stype ) * n );
}

void NoCommunicator::gatherImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    common::scalar::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )

    memcpy( allVals, myVals, common::typeSize( stype ) * n );
}

void NoCommunicator::gatherVImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const IndexType sizes[],
    common::scalar::ScalarType stype ) const
{
    SCAI_ASSERT_EQ_ERROR( root, 0, "" )
    SCAI_ASSERT_EQ_ERROR( sizes[0], n, "" )

    memcpy( allVals, myVals, common::typeSize( stype ) * n );
}

void NoCommunicator::swapImpl( void*, const IndexType, const PartitionId partner, const common::scalar::ScalarType ) const
{
    SCAI_ASSERT_EQ_ERROR( partner, 0, "" )
}

/* --------------------------------------------------------------- */

void NoCommunicator::getProcessorName( char* name ) const
{
    size_t len = maxProcessorName();

    memset( name, 0, len * sizeof( char ) );

    gethostname( name,len );
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

void NoCommunicator::writeAt( std::ostream& stream ) const
{
    stream << "NoComm";
}

/* --------------------------------------------------------------- */

static common::weak_ptr<class NoCommunicator> theNoCommunicatorInstance;

CommunicatorPtr NoCommunicator::create()
{
    common::shared_ptr<NoCommunicator> communicator;

    // use the last communicatorInstance if it is still valid

    if ( theNoCommunicatorInstance.expired() )
    {
        // create a new instance of NoCommunicator and keep it for further uses
        communicator = common::shared_ptr<NoCommunicator>( new NoCommunicator() );
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

Communicator::CommunicatorKind NoCommunicator::createValue()
{
    return NO;
}

} /* end namespace dmemo */

} /* end namespace scai */
