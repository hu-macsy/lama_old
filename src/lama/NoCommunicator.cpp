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
    : Communicator( "none" )
{
    LAMA_LOG_DEBUG( logger, "NoCommunicator()" )
}

NoCommunicator::~NoCommunicator()
{
    LAMA_LOG_DEBUG( logger, "~NoCommunicator()" )
}

ContextPtr NoCommunicator::getCommunicationContext() const
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

void NoCommunicator::all2all( int* recvValues, const int* sendValues ) const
{
    recvValues[0] = sendValues[0];
}

void NoCommunicator::exchangeByPlanImpl(
    void* const recvData,
    const CommunicationPlan& recvPlan,
    const void* const sendData,
    const CommunicationPlan& sendPlan,
    int elemSize ) const
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

    memcpy( recvData, sendData, quantity * elemSize );
}

void NoCommunicator::exchangeByPlan(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(int) );
}

void NoCommunicator::exchangeByPlan(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(float) );
}

void NoCommunicator::exchangeByPlan(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(double) );
}

SyncToken* NoCommunicator::exchangeByPlanAsync(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(int) );
    return new NoSyncToken();
}

SyncToken* NoCommunicator::exchangeByPlanAsync(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(float) );
    return new NoSyncToken();
}

SyncToken* NoCommunicator::exchangeByPlanAsync(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan, sizeof(double) );
    return new NoSyncToken();
}

IndexType NoCommunicator::shiftImpl(
    double targetVals[],
    const IndexType targetSize,
    const double sourceVals[],
    const IndexType sourceSize,
    const int direction ) const
{
    if ( direction != 0 )
    {
        LAMA_LOG_WARN( logger, "shift<double> for NoCommunicator, dir = " << direction )
    }

    return Communicator::shift0( targetVals, targetSize, sourceVals, sourceSize );
}

IndexType NoCommunicator::shiftImpl(
    float targetVals[],
    const IndexType targetSize,
    const float sourceVals[],
    const IndexType sourceSize,
    const int direction ) const
{
    if ( direction != 0 )
    {
        LAMA_LOG_WARN( logger, "shift<float> for NoCommunicator, dir = " << direction )
    }
    return Communicator::shift0( targetVals, targetSize, sourceVals, sourceSize );
}

IndexType NoCommunicator::shiftImpl(
    int targetVals[],
    const IndexType targetSize,
    const int sourceVals[],
    const IndexType sourceSize,
    const int direction ) const
{
    if ( direction != 0 )
    {
        LAMA_LOG_WARN( logger, "shift<int> for NoCommunicator, dir = " << direction )
    }
    return Communicator::shift0( targetVals, targetSize, sourceVals, sourceSize );
}

void NoCommunicator::maxloc( double&, int&, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::maxloc( float&, int&, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::maxloc( int&, int&, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::bcast( double[], const IndexType, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::bcast( float[], const IndexType, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::bcast( int[], const IndexType, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::bcast( std::string&, const PartitionId root ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
}

void NoCommunicator::scatter( double myvals[], const IndexType n, const PartitionId root, const double allvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::scatter(
    double myvals[],
    const IndexType n,
    const PartitionId root,
    const double allvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::scatter(
    float myvals[],
    const IndexType n,
    const PartitionId root,
    const float allvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::scatter(
    int myvals[],
    const IndexType n,
    const PartitionId root,
    const int allvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = allvals[i];
    }
}

void NoCommunicator::gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::gather(
    double allvals[],
    const IndexType n,
    const PartitionId root,
    const double myvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_EQUAL_ERROR( sizes[0], n )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::gather(
    float allvals[],
    const IndexType n,
    const PartitionId root,
    const float myvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_ERROR( sizes[0] == n, "illegal array sizes, sizes[0] = " << sizes[0] << ", expected = " << n )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::gather(
    int allvals[],
    const IndexType n,
    const PartitionId root,
    const int myvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_EQUAL_ERROR( root, 0 )
    LAMA_ASSERT_ERROR( sizes[0] == n, "illegal array sizes, sizes[0] = " << sizes[0] << ", expected = " << n )

    for ( int i = 0; i < n; i++ )
    {
        allvals[i] = myvals[i];
    }
}

void NoCommunicator::swap( double[], const IndexType, const PartitionId partner ) const
{
    LAMA_ASSERT_EQUAL_ERROR( partner, 0 )
}

void NoCommunicator::swap( float[], const IndexType, const PartitionId partner ) const
{
    LAMA_ASSERT_EQUAL_ERROR( partner, 0 )
}

void NoCommunicator::swap( int[], const IndexType, const PartitionId partner ) const
{
    LAMA_ASSERT_EQUAL_ERROR( partner, 0 )
}

float NoCommunicator::sum( const float value ) const
{
    return value;
}

double NoCommunicator::sum( const double value ) const
{
    return value;
}

int NoCommunicator::sum( const int value ) const
{
    return value;
}

size_t NoCommunicator::sum( const size_t value ) const
{
    return value;
}

float NoCommunicator::min( const float value ) const
{
    return value;
}

float NoCommunicator::max( const float value ) const
{
    return value;
}

double NoCommunicator::min( const double value ) const
{
    return value;
}

double NoCommunicator::max( const double value ) const
{
    return value;
}

int NoCommunicator::min( const int value ) const
{
    return value;
}

int NoCommunicator::max( const int value ) const
{
    return value;
}

void NoCommunicator::gather( vector<IndexType>& values, IndexType value ) const
{
    // build a vector of just a single value

    values.clear();
    values.push_back( value );
}

void NoCommunicator::gather( vector<float>& values, float value ) const
{
    // build a vector of just a single value

    values.clear();
    values.push_back( value );
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

