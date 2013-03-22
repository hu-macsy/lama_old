/**
 * @file OpenShMemInterface.cpp
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
 * @brief OpenShMemInterface.cpp
 * @author Michael Drost
 * @date 03.02.2012
 * $Id$
 */

#include <lama/pgas/OpenShMemInterface.hpp>

#include <mpp/shmem.h>

#ifdef OSH_VERBOSE_DEBUG
#include <iostream>
#include <sstream>
#include <string>
std::ostringstream OpenShMemInterface::mProtocol;
#pragma message("OpenShMemInterface compiled with OSH_VERBOSE_DEBUG set, this has a VERY HUGE performance impact!!!")
#endif
namespace lama
{
LAMA_LOG_DEF_LOGGER( OpenShMemInterface::logger, "PGASInterface.OpenShMemInterface" )
#ifdef LAMA_LOG_TRACE
long OpenShMemInterface::mAllocated = 0;
long OpenShMemInterface::mBarrierNum = 0;
#endif
OpenShMemInterface::OpenShMemInterface()
{
    start_pes( 0 );
    mRank = _my_pe();
    mSize = _num_pes();
    pSync = static_cast<long*>( OpenShMemInterface::allocate( _SHMEM_BCAST_SYNC_SIZE * sizeof(long) ) );
    pWork = static_cast<void*>( OpenShMemInterface::allocate( _SHMEM_REDUCE_SYNC_SIZE * sizeof(double) ) );
}

OpenShMemInterface::~OpenShMemInterface()
{
    OpenShMemInterface::free( pSync, _SHMEM_BCAST_SYNC_SIZE * sizeof(long) );
    OpenShMemInterface::free( pWork, _SHMEM_REDUCE_SYNC_SIZE * sizeof(double) );
#ifdef LAMA_LOG_TRACE
    if ( mAllocated != 0 )
    {
        LAMA_LOG_WARN( logger, "There are still " << mAllocated << " B allocated in the OpenShMem pinned memory!" )
    }
#endif
}

PartitionId OpenShMemInterface::getRank() const
{
    return mRank;
}

PartitionId OpenShMemInterface::getSize() const
{
    return mSize;
}

void *OpenShMemInterface::allocate( size_t size ) const
{
#ifdef LAMA_LOG_TRACE
    OpenShMemInterface::mAllocated += size * 2;
#endif
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Will allocate " << size*2 << "B currently allocated: " << mAllocated )
    syncronizeAll();
    void* ptr = shmalloc( size * 2 );
    syncronizeAll();
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Allocated " << size*2 << "B @ " << ptr << " currently allocated: " << mAllocated )
#ifdef OSH_VERBOSE_DEBUG
    mProtocol << "Allocated " << size*2 << "B @ " << ptr << " currently allocated: " << mAllocated << std::endl;
    Exception::addCallStack(mProtocol);
    if(ptr==NULL)
    {
        std::cout << "#############################################################################" << std::endl;
        std::cout << "##                          BEGIN ALLOCATIONS DUMP                         ##" << std::endl;
        std::cout << "#############################################################################" << std::endl;
        std::cout << std::endl << mProtocol.str() << std::endl << std::endl;
        std::cout << "#############################################################################" << std::endl;
        std::cout << "##                          END ALLOCATIONS DUMP                           ##" << std::endl;
        std::cout << "#############################################################################" << std::endl;
    }
#endif
    LAMA_ASSERT( ptr != NULL, "Openshmem could not allocate on Rank " << getRank() )
    return ptr;
}

void OpenShMemInterface::free( void *ptr, const size_t size ) const
{
#ifdef LAMA_LOG_TRACE
    OpenShMemInterface::mAllocated -= size * 2;
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Freed " << size << "B @ " << ptr << " currently allocated: " << mAllocated )
//    LAMA_ASSERT(mAllocated!=1588,"Here it is ... ");
#endif
#ifdef OSH_VERBOSE_DEBUG
    mProtocol << "Freed " << size*2 << "B @ " << ptr << " currently allocated: " << mAllocated << std::endl;
    Exception::addCallStack(mProtocol);
#endif
    shfree( ptr );
    return;
}

void OpenShMemInterface::syncronizeAll() const
{
#ifdef LAMA_LOG_TRACE
    mBarrierNum++;
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Barrier " << mBarrierNum << " currently allocated: " << mAllocated )
//    LAMA_ASSERT(mBarrierNum!=20,"Here it is ... ");
#endif
#ifdef OSH_VERBOSE_DEBUG
    mProtocol << "Barrier " << mBarrierNum << " currently allocated: " << mAllocated << std::endl;
    Exception::addCallStack(mProtocol);
#endif
    shmem_barrier_all();
}

bool OpenShMemInterface::isPinned( const void * const ptr ) const
{
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": isPinned " << ptr << " (" << shmem_addr_accessible(const_cast<void*>(ptr), 0) << ")" )
    return shmem_addr_accessible( const_cast<void*>( ptr ), 0 )
}

PGASSyncToken* OpenShMemInterface::getSyncToken( int ) const
{
    return new PGASSyncToken();
}

PGASCommunicationKind OpenShMemInterface::getPreferredCommunicationKind() const
{
    return PGASget;
}

void OpenShMemInterface::get( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Fetching " << length << "B from " << srcpe << ".  currently allocated: " << mAllocated )
    shmem_getmem( static_cast<void*>( dst ), static_cast<const void*>( src ), length, srcpe )
}

void OpenShMemInterface::put( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_LOG_TRACE( logger,
                    "Rank: " << getRank() << ": Putting " << length << "B to " << srcpe << ".  currently allocated: " << mAllocated )
    shmem_putmem( static_cast<void*>( dst ), static_cast<const void*>( src ), length, srcpe )
}

void OpenShMemInterface::writeAt( std::ostream & stream ) const
{
    stream << "OpenShMemInterface Rank: " << getRank() << "/" << getSize();
}

void OpenShMemInterface::initSync() const
{
    for ( int i = 0; i < _SHMEM_BCAST_SYNC_SIZE; ++i )
    {
        pSync[i] = _SHMEM_SYNC_VALUE;
    }
}

}
