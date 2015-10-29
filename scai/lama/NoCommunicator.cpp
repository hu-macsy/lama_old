/**
 * @file NoCommunicator.cpp
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
 * @brief NoCommunicator.cpp
 * @author Thomas Brandes
 * @date 15.03.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/NoCommunicator.hpp>

// local library
#include <scai/lama/CommunicationPlan.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/Assert.hpp>
#include <scai/common/weak_ptr.hpp>

using namespace std;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( NoCommunicator::logger, "Communicator.NoCommunicator" )

NoCommunicator::NoCommunicator()
                : CRTPCommunicator<NoCommunicator>( "none" )
{
    SCAI_LOG_DEBUG( logger, "NoCommunicator()" )
}

NoCommunicator::~NoCommunicator()
{
    SCAI_LOG_DEBUG( logger, "~NoCommunicator()" )
}

hmemo::ContextPtr NoCommunicator::getCommunicationContext( const hmemo::ContextArray& ) const
{
    return hmemo::Context::getContextPtr( hmemo::context::Host );
}

bool NoCommunicator::isEqual( const Communicator& other ) const
{
    return typeid( *this ) == typeid( other );
}

Communicator::ThreadSafetyLevel NoCommunicator::getThreadSafetyLevel() const
{
    return Communicator::Multiple;
}

PartitionId NoCommunicator::getSize() const
{
    return 1;
}

PartitionId NoCommunicator::getRank() const
{
    return 0;
}

PartitionId NoCommunicator::getNodeSize() const
{
    return 1;
}

PartitionId NoCommunicator::getNodeRank() const
{
    return 0;
}

void NoCommunicator::all2all( IndexType recvValues[], const IndexType sendValues[] ) const
{
    recvValues[0] = sendValues[0];
}

/* ---------------------------------------------------------------------------------- */
/*      exchangeByPlan                                                                */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void NoCommunicator::exchangeByPlanImpl(
    ValueType recvData[],
    const CommunicationPlan& recvPlan,
    const ValueType sendData[],
    const CommunicationPlan& sendPlan ) const
{
    SCAI_ASSERT_EQUAL_ERROR( recvPlan.size(), sendPlan.size() )

    if( 0 == recvPlan.size() && 0 == sendPlan.size() )
    {
        return;
    }

    // send / recv plan have maximal one value

    SCAI_ASSERT_EQUAL_ERROR( 1, recvPlan.size() )

    int quantity = recvPlan[0].quantity;

    // recv and send plan must have same quantity

    SCAI_ASSERT_EQUAL_ERROR( quantity, sendPlan[0].quantity )

    // self copy of send data to recv data

    memcpy( recvData, sendData, quantity * sizeof(ValueType) );
}

template<typename ValueType>
tasking::SyncToken* NoCommunicator::exchangeByPlanAsyncImpl(
    ValueType recvData[],
    const CommunicationPlan& recvPlan,
    const ValueType sendData[],
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
    return new tasking::NoSyncToken();
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
IndexType NoCommunicator::shiftImpl(
    ValueType[],
    const IndexType,
    const PartitionId,
    const ValueType[],
    const IndexType,
    const PartitionId ) const
{
    COMMON_THROWEXCEPTION( "shiftImpl should never be called for NoCommunicator" )
}

template<typename ValueType>
tasking::SyncToken* NoCommunicator::shiftAsyncImpl(
    ValueType[],
    const PartitionId,
    const ValueType[],
    const PartitionId,
    const IndexType ) const
{
    COMMON_THROWEXCEPTION( "shiftAsyncImpl should never be called for NoCommunicator" )
}

template<typename ValueType>
void NoCommunicator::maxlocImpl( ValueType&, IndexType&, const PartitionId root ) const
{
    // nothing to do
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )
}

template<typename ValueType>
void NoCommunicator::bcastImpl( ValueType[], const IndexType, const PartitionId root ) const
{
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )
}

template<typename ValueType>
void NoCommunicator::all2allvImpl( ValueType**, IndexType[], ValueType**, IndexType[] /*sendCount[]*/ ) const
{
   // SCAI_ASSERT_EQUAL_ERROR(sendCount,0)
}

template<typename ValueType>
void NoCommunicator::scatterImpl(
    ValueType myvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType allvals[] ) const
{
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )

    for( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

template<typename ValueType>
void NoCommunicator::scatterVImpl(
    ValueType myvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType allvals[],
    const IndexType sizes[] ) const
{
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )
    SCAI_ASSERT_EQUAL_ERROR( sizes[0], n )

    for( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

template<typename ValueType>
void NoCommunicator::gatherImpl(
    ValueType allvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType myvals[] ) const
{
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )

    for( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

template<typename ValueType>
void NoCommunicator::gatherVImpl(
    ValueType allvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType myvals[],
    const IndexType sizes[] ) const
{
    SCAI_ASSERT_EQUAL_ERROR( root, 0 )
    SCAI_ASSERT_EQUAL_ERROR( sizes[0], n )

    for( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

template<typename ValueType>
void NoCommunicator::swapImpl( ValueType[], const IndexType, const PartitionId partner ) const
{
    SCAI_ASSERT_EQUAL_ERROR( partner, 0 )
}

template<typename ValueType>
ValueType NoCommunicator::sumImpl( const ValueType value ) const
{
    return value;
}

template<typename ValueType>
ValueType NoCommunicator::minImpl( const ValueType value ) const
{
    return value;
}

template<typename ValueType>
ValueType NoCommunicator::maxImpl( const ValueType value ) const
{
    return value;
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

    if( theNoCommunicatorInstance.expired() )
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

std::string NoCommunicator::createValue()
{
    return "none";
}

} /* end namespace lama */

} /* end namespace scai */
