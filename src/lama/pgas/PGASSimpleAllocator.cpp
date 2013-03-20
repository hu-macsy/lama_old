/**
 * @file PGASSimpleAllocator.cpp
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
 * @brief PGASSimpleAllocator.cpp
 * @author Michael Drost
 * @date 06.03.2012
 * $Id$
 */

#include <lama/pgas/PGASSimpleAllocator.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <iostream>

namespace lama
{

PGASSimpleAllocator::PGASSimpleAllocator( void* basepointer, size_t size )
{
    mBasePointer = basepointer;
    mActualPointer = mBasePointer;
    mSize = size;
    mFreeSpace = mSize;
    //mLastPtrs.push_back(mBasePointer);
//    std::cout << "SPACE LEFT:" << mFreeSpace << std::endl;
}

PGASSimpleAllocator::~PGASSimpleAllocator()
{
    // TODO Auto-generated destructor stub
}

void* PGASSimpleAllocator::allocate( size_t size )
{
    char* temp = static_cast<char*>( mActualPointer );
    mActualPointer = temp + size;
//    std::cout << "SPACE LEFT:" << mFreeSpace << std::endl;
//    std::cout << "SIZE:" << size << std::endl;
    mFreeSpace -= size;
//    std::cout << "SPACE LEFT:" << mFreeSpace << std::endl;
//    std::cout << "Allocated @ " << static_cast<void*>(temp) << std::endl;
    LAMA_ASSERT( mFreeSpace > 0, "Out of Memory" );
    mLastPtrs.push_back( mActualPointer );
    return temp;
}

size_t PGASSimpleAllocator::getOffset( const void* actualptr )
{
    const char* base = static_cast<const char*>( mBasePointer );
    const char* ptr = static_cast<const char*>( actualptr );
    return ptr - base;
}

void PGASSimpleAllocator::recalcFreeSpace()
{
    if( mLastPtrs.size() == 0 )
    {
        mFreeSpace = mSize;
        return;
    }
    mFreeSpace = mSize - ( static_cast<char*>( mActualPointer ) - static_cast<char*>( mBasePointer ) );
}

void PGASSimpleAllocator::free( void* ptr, size_t size )
{
    //Nothing to do here
//    std::cout << "freeing SPACE LEFT:" << mFreeSpace << std::endl;
    char* temp = static_cast<char*>( ptr ) + size;
    if( ptr != mBasePointer )
    {
        mLastPtrs.remove( temp );
    }
    mActualPointer = mLastPtrs.back();
    recalcFreeSpace();
//    std::cout << "freed SPACE LEFT:" << mFreeSpace << std::endl;
}

bool PGASSimpleAllocator::isAllocated( const void* ptr )
{
    return ( ptr >= mBasePointer && ptr < mActualPointer );
}
}
