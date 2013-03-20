/**
 * @file PGASSimpleAllocator.hpp
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
 * @brief PGASSimpleAllocator.hpp
 * @author Michael Drost
 * @date 06.03.2012
 * $Id$
 */
#ifndef LAMA_PGASSIMPLEALLOCATOR_HPP_
#define LAMA_PGASSIMPLEALLOCATOR_HPP_

#include <lama/pgas/PGASAllocator.hpp>

#include <list>

namespace lama
{

class PGASSimpleAllocator: public lama::PGASAllocator
{
private:
    void* mBasePointer;
    size_t mSize;
    void* mActualPointer;
    std::list<void*> mLastPtrs;
    size_t mFreeSpace;
public:
    PGASSimpleAllocator( void* basepointer, size_t size );
    virtual ~PGASSimpleAllocator();
    virtual void* allocate( size_t size );
    virtual void free( void* ptr, size_t size );
    void recalcFreeSpace();
    virtual bool isAllocated( const void* ptr );
    virtual size_t getOffset( const void* ptr );
};

} /* namespace lama */

#endif // LAMA_PGASSIMPLEALLOCATOR_HPP_
