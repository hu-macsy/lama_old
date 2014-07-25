/**
 * @file NoCommunicator.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#include <lama/NoCommunicator.hpp>

// others
#include <lama/NoSyncToken.hpp>
#include <lama/CommunicationPlan.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

using namespace std;

namespace lama
{

LAMA_LOG_DEF_LOGGER( NoCommunicator::logger, "Communicator.NoCommunicator" )

NoCommunicator::NoCommunicator()
    : CRTPCommunicator<NoCommunicator>( "none" )
{
    LAMA_LOG_DEBUG( logger, "NoCommunicator()" )
}

NoCommunicator::~NoCommunicator()
{
    LAMA_LOG_DEBUG( logger, "~NoCommunicator()" )
}

ContextPtr NoCommunicator::getCommunicationContext( const _LAMAArray& ) const
{
    return ContextFactory::getContext( Context::Host );
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

template<typename T>
void NoCommunicator::exchangeByPlanImpl(
    T recvData[],
    const CommunicationPlan& recvPlan,
    const T sendData[],
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_EQUAL_ERROR( recvPlan.size(), sendPlan.size() )

    if ( 0 == recvPlan.size() && 0 == sendPlan.size() )
    {
        return;
    }

    // send / recv plan have maximal one value 

    LAMA_ASSERT_EQUAL_ERROR( 1, recvPlan.size() )

    int quantity = recvPlan[0].quantity;

    // recv and send plan must have same quantity

    LAMA_ASSERT_EQUAL_ERROR( quantity, sendPlan[0].quantity )

    // self copy of send data to recv data

    memcpy( recvData, sendData, quantity * sizeof( T ) );
}

template<typename T>
SyncToken* NoCommunicator::exchangeByPlanAsyncImpl(
    T recvData[],
    const CommunicationPlan& recvPlan,
    const T sendData[],
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
    return new NoSyncToken();
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename T>
IndexType NoCommunicator::shiftImpl(
    T[],
    const IndexType,
    const PartitionId,
    const T[],
    const IndexType,
    const PartitionId ) const
{
    LAMA_THROWEXCEPTION( "shiftImpl should never be called for NoCommunicator" )
}

template<typename T>
SyncToken* NoCommunicator::shiftAsyncImpl(
    T[],
    const PartitionId,
    const T[],
    const PartitionId,
    const IndexType ) const
{
    LAMA_THROWEXCEPTION( "shiftAsyncImpl should never be called for NoCommunicator" )
}

template<typename T>
void NoCommunicator::maxlocImpl( T&, IndexType&, const PartitionId root ) const
{
    // nothing to do
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}
 
template<typename T>
void NoCommunicator::bcastImpl( T[], const IndexType, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

template<typename T>
void NoCommunicator::scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

template<typename T>
void NoCommunicator::scatterVImpl(
    T myvals[],
    const IndexType n,
    const PartitionId root,
    const T allvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

template<typename T>
void NoCommunicator::gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

template<typename T>
void NoCommunicator::gatherVImpl(
    T allvals[],
    const IndexType n,
    const PartitionId root,
    const T myvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

template<typename T>
void NoCommunicator::swapImpl( T[], const IndexType, const PartitionId partner ) const
{
    LAMA_ASSERT_EQUAL_ERROR( partner, 0 )
}

template<typename T>
T NoCommunicator::sumImpl( const T value ) const
{
    return value;
}

template<typename T>
T NoCommunicator::minImpl( const T value ) const
{
    return value;
}

template<typename T>
T NoCommunicator::maxImpl( const T value ) const
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

}

