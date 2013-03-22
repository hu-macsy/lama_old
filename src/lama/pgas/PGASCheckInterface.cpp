/**
 * @file PGASCheckInterface.cpp
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
 * @brief PGASCheckInterface.cpp
 * @author Michael Drost
 * @date 07.02.2012
 * $Id$
 */

#include <lama/pgas/PGASCheckInterface.hpp>

#include <lama/exception/LAMAAssert.hpp>

namespace lama
{
void *PGASCheckInterface::allocate( size_t size ) const
{
    return mInterface->allocate( size );
}

void PGASCheckInterface::free( void *ptr, const size_t size ) const
{
    LAMA_ASSERT( mInterface->isPinned(ptr), "Memory is not pinned" )
    mInterface->free( ptr, size )
}

void PGASCheckInterface::syncronizeAll() const
{
    mInterface->syncronizeAll();
}

bool PGASCheckInterface::isPinned( const void * const ptr ) const
{
    return mInterface->isPinned( ptr );
}

PGASSyncToken *PGASCheckInterface::getSyncToken( int arg1 ) const
{
    return mInterface->getSyncToken( arg1 );
}

PGASCommunicationKind PGASCheckInterface::getPreferredCommunicationKind() const
{
    return mInterface->getPreferredCommunicationKind();
}

void PGASCheckInterface::get( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    mInterface->get( dst, src, length, srcpe )

}

void PGASCheckInterface::put( void *dst, const void *src, size_t length, int srcpe ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    mInterface->put( dst, src, length, srcpe )
}

std::auto_ptr<SyncToken> lama::PGASCheckInterface::getAsync(
    void *dst,
    const void *src,
    size_t length,
    int srcPE ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    return mInterface->getAsync( dst, src, length, srcPE )
}

std::auto_ptr<SyncToken> lama::PGASCheckInterface::putAsync(
    void *dst,
    const void *src,
    size_t length,
    int srcPE ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    return mInterface->putAsync( dst, src, length, srcPE )
}

std::auto_ptr<SyncToken> lama::PGASCheckInterface::shift(
    void *dst,
    const void *src,
    size_t size,
    PartitionId destRank,
    PartitionId srcRank ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    return mInterface->shift( dst, src, size, destRank, srcRank )
}

std::auto_ptr<SyncToken> lama::PGASCheckInterface::broadcast(
    void *dst,
    const void *src,
    size_t length,
    int srcPE ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    return mInterface->broadcast( dst, src, length, srcPE )
}

std::auto_ptr<SyncToken> lama::PGASCheckInterface::all2all( void *dst, const void *src, size_t elemSize ) const
{
    LAMA_ASSERT( mInterface->isPinned(dst), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(src), "Memory is not pinned" )
    return mInterface->all2all( dst, src, elemSize )
}

void lama::PGASCheckInterface::swap( void *val, const size_t n, const PartitionId partner ) const
{
    LAMA_ASSERT( mInterface->isPinned(val), "Memory is not pinned" )
    mInterface->swap( val, n, partner )
}

void lama::PGASCheckInterface::scatter(
    void *myvals,
    const size_t partSize,
    const PartitionId root,
    const void *allvals ) const
{
    LAMA_ASSERT( mInterface->isPinned(myvals), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(allvals), "Memory is not pinned" )
    mInterface->scatter( myvals, partSize, root, allvals )
}

void lama::PGASCheckInterface::scatter(
    void *myvals,
    const size_t elemSize,
    const PartitionId root,
    const void *allvals,
    const IndexType sizes[] ) const
{
    LAMA_ASSERT( mInterface->isPinned(myvals), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(allvals), "Memory is not pinned" )
    mInterface->scatter( myvals, elemSize, root, allvals, sizes )
}

void lama::PGASCheckInterface::gather(
    void *allvals,
    const size_t partSize,
    const PartitionId root,
    const void *myvals ) const
{
    LAMA_ASSERT( mInterface->isPinned(myvals), "Memory is not pinned" )
    LAMA_ASSERT( mInterface->isPinned(allvals), "Memory is not pinned" )
    mInterface->gather( allvals, partSize, root, myvals )
}

void lama::PGASCheckInterface::gather(
    void *allvals,
    const size_t elemSize,
    const PartitionId root,
    const void *myvals,
    const IndexType sizes[] ) const
{
    LAMA_ASSERT( mInterface->isPinned(myvals), "Memory is not pinned: " << myvals )
    LAMA_ASSERT( mInterface->isPinned(allvals), "Memory is not pinned: " << allvals )
    mInterface->gather( allvals, elemSize, root, myvals, sizes )
}

lama::PGASCheckInterface::PGASCheckInterface( PGASInterface* interface )
    : mInterface( interface )
{
}

PartitionId PGASCheckInterface::getRank() const
{
    return mInterface->getRank();
}

PartitionId PGASCheckInterface::getSize() const
{
    return mInterface->getSize();
}

void PGASCheckInterface::writeAt( std::ostream & stream ) const
{
    stream << "PGASCheckInterface: " << *mInterface.get();
}

PGASCheckInterface::~PGASCheckInterface()
{
}

}

