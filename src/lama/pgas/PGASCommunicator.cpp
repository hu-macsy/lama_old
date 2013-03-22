/**
 * @file PGASCommunicator.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief PGASCommunicator.cpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASCommunicator.hpp>
#include <lama/pgas/GPIInterface.hpp>

namespace lama
{

const int PGASCommunicator::defaultTag = 1;

LAMA_LOG_DEF_LOGGER( PGASCommunicator::logger, "Communicator.PGASCommunicator" )

PGASCommunicator::PGASCommunicator( int&, char** & )
    : Communicator( "PGAS" ), mPGASInterface( PGASInterface::getInstance() ), mThreadSafetyLevel(
        Communicator::Funneled )
{
}

PGASCommunicator::~PGASCommunicator()
{

}

bool PGASCommunicator::isEqual( const Communicator& other ) const
{
    return ( NULL != dynamic_cast<const PGASCommunicator*>( &other ) );
}

Communicator::ThreadSafetyLevel PGASCommunicator::getThreadSafetyLevel() const
{
    return mThreadSafetyLevel;
}

PartitionId PGASCommunicator::getSize() const
{
    return mPGASInterface->getSize();
}

PartitionId PGASCommunicator::getRank() const
{
    return mPGASInterface->getRank();
}

/* ---------------------------------------------------------------------------------- */

void PGASCommunicator::all2all( int* recvSizes, const int* sendSizes ) const
{
    mPGASInterface->all2all( static_cast<void*>( recvSizes ), static_cast<const void*>( sendSizes ), sizeof(int) );
}

/* ---------------------------------------------------------------------------------- */

void PGASCommunicator::exchangeByPlan(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

void PGASCommunicator::exchangeByPlan(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

void PGASCommunicator::exchangeByPlan(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> PGASCommunicator::exchangeByPlanAsync(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> PGASCommunicator::exchangeByPlanAsync(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> PGASCommunicator::exchangeByPlanAsync(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

template<typename T>
void PGASCommunicator::exchangeByPlanImpl(
    T* const recvData,
    const CommunicationPlan& recvPlan,
    const T* const sendData,
    const CommunicationPlan& sendPlan ) const
{

    exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
//    LAMA_REGION( "PGASCommunicator::exchangeByPlan" )
//    LAMA_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
//    LAMA_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
//    LAMA_LOG_INFO( logger, *this << ": exchange for values of type " << typeid( T ).name()
//              << ", send to " << sendPlan.size()
//              << " processors, recv from " << recvPlan.size() )
//    int maxReceives = recvPlan.size();
//    int noReceives  = 0;  // will be incremented
//    T* recvDataForMe = NULL;
//    IndexType recvDataForMeSize = 0;
//    boost::scoped_array<PGAS_Request> commRequest( new PGAS_Request[maxReceives] );
//
//    // setup receives for each entry in receive plan
//
//    for ( PartitionId i = 0; i < maxReceives; ++i )
//    {
//        IndexType quantity = recvPlan[i].quantity;
//        T* recvDataForI = recvData + recvPlan[i].offset;
//        PartitionId p = recvPlan[i].partitionId;
//        LAMA_LOG_DEBUG( logger, *this << ": receive " << quantity << " elements"
//                   << " from processor " << p )
//
//        if ( p != mRank )
//        {
//            commRequest[noReceives] = startrecv( recvDataForI, quantity, p );
//            noReceives++;
//        }
//        else
//        {
//            recvDataForMe = recvDataForI;
//            recvDataForMeSize = quantity;
//        }
//    }
//
//    // send the data via sendPlan
//
//    for ( PartitionId i = 0; i < sendPlan.size(); ++i )
//    {
//        IndexType quantity = sendPlan[i].quantity;
//        const T* sendDataForI = sendData + sendPlan[i].offset;
//        PartitionId p = sendPlan[i].partitionId;
//        LAMA_LOG_DEBUG( logger, *this << ": send " << quantity << " elements"
//                   << " to processor " << p )
//
//        if ( p != mRank )
//        {
//            send( sendDataForI, quantity, p );
//        }
//        else
//        {
//            LAMA_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
//            LAMA_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )
//
//            for ( IndexType k = 0; k < recvDataForMeSize; k++ )
//            {
//                recvDataForMe[k] = sendDataForI[k];
//            }
//        }
//    }
//
//    // wait for completion of receives
//    boost::scoped_array<PGAS_Status> statuses( new PGAS_Status[noReceives] );
//    PGASCALL( logger, PGAS_Waitall( noReceives, commRequest.get(), statuses.get() ), "PGAS_Waitall" );
//    // ToDo: check for correct sizes, was done in earlier version, but is now redundant
}

/* ---------------------------------------------------------------------------------- */

template<typename T>
std::auto_ptr<SyncToken> PGASCommunicator::exchangeByPlanAsyncImpl(
    T* const recvData,
    const CommunicationPlan& recvPlan,
    const T* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    LAMA_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    LAMA_LOG_INFO( logger,
                   *this << ": exchange for values of type " << typeid( T ).name() << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )

    //TODO Optimize
    //While Using PGAS for exchangeByPlan we need the offsets of our Data on the
    //remote nodes. so we need to communicate with every node to get our offsets on that node.
    int* offsets = static_cast<int*>( mPGASInterface->allocate( getSize() * sizeof(int) ) );
    {
        //Temporary array for the local offsets copied from our sendPlan
        int* temp = static_cast<int*>( mPGASInterface->allocate( getSize() * sizeof(int) ) );
        //this temporary Array could stay uninitialized because only nodes we have values for
        //are accessing the values.
//        for(int i = 0; i<getSize(); ++i)
//        {
//            temp[i] = -1;
//        }
        //Copy our offsets from the sendPlan
        for ( int i = 0; i < sendPlan.size(); ++i )
        {
            temp[sendPlan[i].partitionId] = sendPlan[i].offset;
        }

        //Communicate with the other nodes to get the offsets.
        all2all( offsets, temp );
        //Free the temporary array
        mPGASInterface->free( temp, getSize() * sizeof(int) );
    }
#ifdef LAMA_LOG_TRACE
    for ( int i = 0; i < sendPlan.size(); ++i )
    {
        LAMA_LOG_TRACE( logger,
                        "Sending " << sendPlan[i].quantity*sizeof(T) << " Bytes of data to " << sendPlan[i].partitionId << "from" << sendPlan[i].offset )
    }
#endif

    std::auto_ptr<SyncToken> token( mPGASInterface->getSyncToken( 0 ) )
    for ( int i = 0; i < recvPlan.size(); ++i )
    {
        PartitionId p = recvPlan[i].partitionId;
        const T* sendDataForI = sendData; //+ offsets[p];
        T* recvDataForI = recvData + recvPlan[i].offset;
        IndexType quantity = recvPlan[i].quantity;
        token.reset( ( mPGASInterface->getAsync( recvDataForI, sendDataForI, quantity * sizeof(T), p ) ).release() );
        LAMA_LOG_TRACE( logger, "Reciving " << quantity*sizeof(T) << " Bytes of data from " << p << "@" << offsets[p] )
    }
    //Free our offsets array. This can be done asynchronous because no remote node will access it.
    mPGASInterface->free( offsets, getSize() * sizeof(int) )

    return token;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename T>
IndexType PGASCommunicator::shiftImpl(
    T recvVals[],
    const IndexType recvSize,
    const PartitionId source,
    const T sendVals[],
    const IndexType,
    const PartitionId dest ) const
{
    synchronize();
    shiftAsyncImpl( recvVals, source, sendVals, dest, recvSize );
    return recvSize;
}

IndexType PGASCommunicator::shift(
    double recvData[],
    const IndexType recvSize,
    const double sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftImpl( recvData, recvSize, source, sendVals, sendSize, dest );
}

IndexType PGASCommunicator::shift(
    float recvData[],
    const IndexType recvSize,
    const float sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftImpl( recvData, recvSize, source, sendVals, sendSize, dest );
}

IndexType PGASCommunicator::shift(
    int recvData[],
    const IndexType recvSize,
    const int sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftImpl( recvData, recvSize, source, sendVals, sendSize, dest );
}

/* ---------------------------------------------------------------------------------- */
/*              shiftAsync                                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
std::auto_ptr<SyncToken> PGASCommunicator::shiftAsyncImpl(
    T recvVals[],
    const PartitionId source,
    const T sendVals[],
    const PartitionId dest,
    const IndexType size ) const
{
    return mPGASInterface->shift( static_cast<void*>( recvVals ), static_cast<const void*>( sendVals ),
                                  size * sizeof(T), dest, source );
}

std::auto_ptr<SyncToken> PGASCommunicator::shiftAsync(
    double recvVals[],
    const double sendVals[],
    const IndexType size,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger, *this << ": shiftAsync size = " << size << ", direction = " << direction )

    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncImpl( recvVals, source, sendVals, dest, size );
}

std::auto_ptr<SyncToken> PGASCommunicator::shiftAsync(
    float recvVals[],
    const float sendVals[],
    const IndexType size,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncImpl( recvVals, source, sendVals, dest, size );
}

std::auto_ptr<SyncToken> PGASCommunicator::shiftAsync(
    int recvVals[],
    const int sendVals[],
    const IndexType size,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncImpl( recvVals, source, sendVals, dest, size );
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

float PGASCommunicator::sum( const float value ) const
{
    return mPGASInterface->sumToAll( value );
}

double PGASCommunicator::sum( const double value ) const
{
    return mPGASInterface->sumToAll( value );
}

size_t PGASCommunicator::sum( const size_t value ) const
{
    return mPGASInterface->sumToAll( value );
}

int PGASCommunicator::sum( const int value ) const
{
    return mPGASInterface->sumToAll( value );
}

/* ---------------------------------------------------------------------------------- */
/*              min / max reduction                                                   */
/* ---------------------------------------------------------------------------------- */

float PGASCommunicator::min( const float value ) const
{
    return mPGASInterface->minToAll( value );
}

double PGASCommunicator::min( const double value ) const
{
    return mPGASInterface->minToAll( value );
}

int PGASCommunicator::min( const int value ) const
{
    return mPGASInterface->minToAll( value );
}

float PGASCommunicator::max( const float value ) const
{
    return mPGASInterface->maxToAll( value );
}

double PGASCommunicator::max( const double value ) const
{
    return mPGASInterface->maxToAll( value );
}

int PGASCommunicator::max( const int value ) const
{
    return mPGASInterface->maxToAll( value );
}

void PGASCommunicator::gather( std::vector<float>& values, float value ) const
{
    // build a vector of just a single value
    values.clear();
    float* temp = static_cast<float*>( mPGASInterface->allocate( mPGASInterface->getSize() * sizeof(float) ) );
    float* myvals = static_cast<float*>( mPGASInterface->allocate( sizeof(float) ) );
    *myvals = value;

    mPGASInterface->gather( temp, 1, 0, myvals );
    mPGASInterface->broadcast( temp, temp, mPGASInterface->getSize() * sizeof(float), 0 );

    mPGASInterface->free( temp, mPGASInterface->getSize() * sizeof(float) );
    mPGASInterface->free( myvals, sizeof(float) );
    values.resize( getSize(), 0.0 );
}

void PGASCommunicator::synchronize() const
{
    mPGASInterface->syncronizeAll();
}

/* ---------------------------------------------------------------------------------- */
/*      bcast                                                                         */
/* ---------------------------------------------------------------------------------- */

void PGASCommunicator::bcast( int val[], const IndexType n, const PartitionId root ) const
{
    mPGASInterface->broadcast( val, val, n * sizeof(int), root );
}

void PGASCommunicator::bcast( double val[], const IndexType n, const PartitionId root ) const
{
    mPGASInterface->broadcast( val, val, n * sizeof(double), root );
}

void PGASCommunicator::bcast( float val[], const IndexType n, const PartitionId root ) const
{
    mPGASInterface->broadcast( val, val, n * sizeof(float), root );
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals )                                           */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const
{
    mPGASInterface->scatter( myvals, n * sizeof(T), root, allvals );

//    Shmemwrapper::barrier_all();
//    LAMA_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
//    LAMA_LOG_DEBUG( logger, *this << ": scatter of " << n
//               << " elements, root = " << root )
//    Shmemwrapper::shgetmem(myvals, (&allvals[n*mRank]), n, root);
//    Shmemwrapper::barrier_all();
}

void PGASCommunicator::scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

void PGASCommunicator::scatter(
    double myvals[],
    const IndexType n,
    const PartitionId root,
    const double allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

void PGASCommunicator::scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::scatterImpl(
    T myvals[],
    const IndexType,
    const PartitionId root,
    const T allvals[],
    const IndexType sizes[] ) const
{
    mPGASInterface->scatter( myvals, sizeof(T), root, allvals, sizes );
//    Shmemwrapper::barrier_all();
//    LAMA_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
//    LAMA_LOG_DEBUG( logger, *this << ": scatter of " << n
//               << " elements, root = " << root )
//
//    int offset = 0;
//    for(int i = 0; i< mRank; ++i) offset+=sizes[i];
//    Shmemwrapper::shgetmem(myvals, (&allvals[offset]), sizes[mRank], root);
//    Shmemwrapper::barrier_all();
}

void PGASCommunicator::scatter(
    float myvals[],
    const IndexType n,
    const PartitionId root,
    const float allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

void PGASCommunicator::scatter(
    double myvals[],
    const IndexType n,
    const PartitionId root,
    const double allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

void PGASCommunicator::scatter(
    int myvals[],
    const IndexType n,
    const PartitionId root,
    const int allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals )                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const
{

    mPGASInterface->gather( allvals, n * sizeof(T), root, myvals );
//    Shmemwrapper::barrier_all();
//    LAMA_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
//    LAMA_LOG_DEBUG( logger, *this << ": gather of " << n
//               << " elements, root = " << root )
//    if(mRank == root)
//        for(int i = 0; i< mSize; ++i)
//        {
//            Shmemwrapper::shgetmem(&allvals[i*n],myvals,n,i);
//        }
//    Shmemwrapper::barrier_all();
}

void PGASCommunicator::gather(
    double allvals[],
    const IndexType n,
    const PartitionId root,
    const double myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

void PGASCommunicator::gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

void PGASCommunicator::gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals, sizes )                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::gatherImpl(
    T allvals[],
    const IndexType,
    const PartitionId root,
    const T myvals[],
    const IndexType sizes[] ) const
{
    mPGASInterface->gather( static_cast<void*>( allvals ), sizeof(T), root, static_cast<const void*>( myvals ), sizes );
}

void PGASCommunicator::gather(
    float allvals[],
    const IndexType n,
    const PartitionId root,
    const float myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

void PGASCommunicator::gather(
    double allvals[],
    const IndexType n,
    const PartitionId root,
    const double myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

void PGASCommunicator::gather(
    int allvals[],
    const IndexType n,
    const PartitionId root,
    const int myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::maxlocImpl( T& val, int& location, PartitionId root ) const
{
    mPGASInterface->maxloc( val, location, root );
//    int vRank = (mRank+mSize-root)%mSize;
//    ShmemPtr in(sizeof(T));
//    ShmemPtr out(sizeof(T));
//    ShmemPtr maxloc(sizeof(PartitionId));
//
//    *maxloc.getTypePtr<int>() = mRank;
//    *out.getTypePtr<T>() = val;
//
//    Shmemwrapper::barrier_all();
//
//    for(int i = (mSize+1)/2; i>0; i/=2)
//    {
//        if(vRank < i)
//        {
//            Shmemwrapper::shgetmem(in.getTypePtr<T>(),out.getTypePtr<T>(),1,(mRank+i)%mSize);
//            if(*out.getTypePtr<T>() < *in.getTypePtr<T>())
//            {
//                Shmemwrapper::shgetmem(maxloc.getTypePtr<int>(),maxloc.getTypePtr<int>(),1,(mRank+i)%mSize);
//                *out.getTypePtr<T>() = *in.getTypePtr<T>();
//            }
//        }
//        Shmemwrapper::barrier_all();
//    }
//    if(mRank == root)
//    {
//        val = *out.getTypePtr<T>();
//        location = *maxloc.getTypePtr<int>();
//    }

}

void PGASCommunicator::maxloc( float& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

void PGASCommunicator::maxloc( double& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

void PGASCommunicator::maxloc( int& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void PGASCommunicator::swapImpl( T val[], const IndexType n, PartitionId partner ) const
{
    if ( partner == mPGASInterface->getRank() )
    {
        return;
    }
    mPGASInterface->swap( val, n * sizeof(T), partner );

//    ShmemPtr temp(sizeof(T)*n);
//    memcpy(temp.getTypePtr<void>(),val,sizeof(T)*n);
//    Shmemwrapper::barrier_all();
//    Shmemwrapper::shgetmem(val,temp.getTypePtr<T>(),n,partner);
//    Shmemwrapper::barrier_all();
}

void PGASCommunicator::swap( double val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

void PGASCommunicator::swap( float val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

void PGASCommunicator::swap( int val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

ContextPtr PGASCommunicator::getCommunicationContext() const
{
//    return ContextFactory::getContext( Context::Openshmem );
    return ContextFactory::getContext( Context::Host );
}

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void PGASCommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful
    stream << "PGAS(" << mPGASInterface->getRank() << ":" << mPGASInterface->getSize() << ")";
}

}
