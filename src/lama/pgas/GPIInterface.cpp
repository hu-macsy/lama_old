/**
 * @file GPIInterface.cpp
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
 * @brief GPIInterface.cpp
 * @author Michael Drost
 * @date 03.02.2012
 * $Id$
 */

#include <lama/pgas/GPIInterface.hpp>
#include <lama/pgas/PGASSimpleAllocator.hpp>

#include <GPI.h>

#include <lama/LAMAArguments.hpp>
#include <lama/exception/LAMAAssert.hpp>

#include <iostream>
#include <omp.h>

namespace lama
{
LAMA_LOG_DEF_LOGGER( GPIInterface::logger, "PGASInterface.GPIInterface" )
#ifdef LAMA_LOG_TRACE
long GPIInterface::mAllocated = 0;
long GPIInterface::mBarrierNum = 0;
#endif

GPIInterface::GPIInterface()
{
    size_t KiB = 1024;
    size_t MiB = KiB * 1024;
    size_t GiB = MiB * 1024;
    size_t size = 6 * GiB; //TODO get from environment
//    if(checkEnv(0,NULL) != 0)
//        LAMA_LOG_WARN(logger,"GPI checkEnv failed!");
//    char * null = "/bully/GPI/mdrost/lama/tests/lama_unit_test PGASCommunicatorTest";

//#pragma omp single
//    {
    std::cout << "STARTING GPI" << std::endl << std::flush;
    LAMA_ASSERT( !(startGPI(LAMAArguments::getArgc(),LAMAArguments::getArgv(),"",size) != 0), "GPI startGPI failed!" )
    std::cout << "STARTED GPI" << std::endl << std::flush;

    mRank = getRankGPI();
    mSize = getNodeCountGPI();
    LAMA_LOG_WARN( logger, "GPIInterface started: Rank: "<< mRank << "/" << mSize )

    mAllocator = new PGASSimpleAllocator( getDmaMemPtrGPI(), size );
//    }
//    start_pes(0);
//    mRank = _my_pe();
//    mSize = _num_pes();
//    pSync = static_cast<long*>(GPIInterface::allocate(_SHMEM_BCAST_SYNC_SIZE * sizeof (long )));
//    pWork = static_cast<void*>(GPIInterface::allocate(_SHMEM_REDUCE_SYNC_SIZE * sizeof (double)));
}

GPIInterface::~GPIInterface()
{
    delete mAllocator;
    syncronizeAll();
    shutdownGPI();
    std::cout << "DOWNED GPI !!!" << std::endl;
//    GPIInterface::free(pSync, _SHMEM_BCAST_SYNC_SIZE * sizeof (long ));
//    GPIInterface::free(pWork, _SHMEM_REDUCE_SYNC_SIZE * sizeof (double ));
//#ifdef LAMA_LOG_TRACE
//    if(mAllocated != 0)
//        LAMA_LOG_WARN(logger,"There are still " << mAllocated << " B allocated in the OpenShMem pinned memory!");
//#endif
}

PartitionId GPIInterface::getRank() const
{
    return mRank;
}

PartitionId GPIInterface::getSize() const
{
    return mSize;
}

void *GPIInterface::allocate( size_t size ) const
{
//    syncronizeAll();
    void* ptr = mAllocator->allocate( size );
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Allocated " << size << "B @ " << mAllocator->getOffset(ptr) << ". currently allocated " << (mAllocated+=size) )
    LAMA_ASSERT( ptr!=NULL, "GPI could not allocate on Rank " << getRank() )
    return ptr;
}

void GPIInterface::free( void *ptr, const size_t size ) const
{
    mAllocator->free( ptr, size );
//    LAMA_ASSERT(0!=0,"Here it is ... ");
//#ifdef LAMA_LOG_TRACE
//    GPIInterface::mAllocated -= size*2;
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Freed " << size << "B @ " << mAllocator->getOffset(ptr) << ". currently allocated " << (mAllocated-=size) )
////    LAMA_ASSERT(mAllocated!=1588,"Here it is ... ");
//#endif
//#ifdef OSH_VERBOSE_DEBUG
//    mProtocol << "Freed " << size*2 << "B @ " << ptr << " currently allocated: " << mAllocated << std::endl;
//    Exception::addCallStack(mProtocol);
//#endif
//    shfree(ptr);
//    return;
}

void GPIInterface::syncronizeAll() const
{
//    static int num = 0;
    LAMA_LOG_TRACE( logger, "Rank: " << getRank() << ": Barrier " << mBarrierNum++ )
//    double start = omp_get_wtime();
    barrierGPI();
//    double end = omp_get_wtime();
    //LAMA_ASSERT(num!=124, "HERE");
//    std::cout << "waited[" << num ++ << "]: " << (end-start)*1000 << "ms" << std::endl;

    //LAMA_ASSERT((end-start)*1000 < 3, "Waiting very long");
//#ifdef LAMA_LOG_TRACE
//    mBarrierNum++;
//    LAMA_LOG_TRACE(logger,"Rank: " << getRank() << ": Barrier " << mBarrierNum << " currently allocated: " << mAllocated)
////    LAMA_ASSERT(mBarrierNum!=20,"Here it is ... ");
//#endif
//#ifdef OSH_VERBOSE_DEBUG
//    mProtocol << "Barrier " << mBarrierNum << " currently allocated: " << mAllocated << std::endl;
//    Exception::addCallStack(mProtocol);
//#endif
//    shmem_barrier_all();
}

bool GPIInterface::isPinned( const void * const ptr ) const
{
    LAMA_LOG_TRACE( logger, "Rank: " << getRank() << ": isPinned " << ptr << " --> " << mAllocator->isAllocated(ptr) )
//    return shmem_addr_accessible(const_cast<void*>(ptr), 0);
    return mAllocator->isAllocated( ptr );
}

PGASSyncToken* GPIInterface::getSyncToken( int ) const
{
    return new PGASSyncToken();
}

PGASCommunicationKind GPIInterface::getPreferredCommunicationKind() const
{
    return PGASget;
}

void GPIInterface::get( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Fetching " << length << "B from " << srcpe << ".  currently allocated: " << mAllocated )
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Fetching " << mAllocator->getOffset(src) << " --> " << mAllocator->getOffset(dst) )
//    syncronizeAll();
    int stat = readDmaGPI( mAllocator->getOffset( dst ), mAllocator->getOffset( src ), length, srcpe, GPIQueue0 );
    LAMA_ASSERT( stat==0, "readDmaGPI failed" )
    waitDmaGPI (GPIQueue0);
//    syncronizeAll();
}

void GPIInterface::put( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Putting " << length << "B to " << srcpe << ".  currently allocated: " << mAllocated )
//    syncronizeAll();
    int stat = writeDmaGPI( mAllocator->getOffset( src ), mAllocator->getOffset( dst ), length, srcpe, GPIQueue0 );
    LAMA_ASSERT( stat==0, "writeDmaGPI failed" )
//    syncronizeAll();
    waitDmaGPI (GPIQueue0);
}

void GPIInterface::writeAt( std::ostream & stream ) const
{
    stream << "GPIInterface Rank: "; // << getRank() << "/" << getSize();
}

}
