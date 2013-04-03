/**
 * @file GPIInterface.hpp
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
 * @brief GPIInterface.hpp
 * @author Michael Drost
 * @date 03.02.2012
 * $Id$
 */
#ifndef LAMA_GPIINTERFACE_HPP_
#define LAMA_GPIINTERFACE_HPP_
//#define OSH_VERBOSE_DEBUG

#include <lama/pgas/PGASInterface.hpp>
#include <lama/pgas/PGASAllocator.hpp>
#include <lama/pgas/PGASSyncToken.hpp>

namespace lama
{
class GPIInterface: public lama::PGASInterface
{
public:
    GPIInterface();
    virtual ~GPIInterface();
    virtual void* allocate( size_t size ) const;
    virtual void free( void* ptr, const size_t size ) const;
    virtual void syncronizeAll() const;
    virtual bool isPinned( const void* const ptr ) const;
    virtual PGASSyncToken* getSyncToken( int arg1 ) const;
    virtual PartitionId getRank() const;
    virtual PartitionId getSize() const;

    virtual PGASCommunicationKind getPreferredCommunicationKind() const;

    virtual void get( void* dst, const void* src, size_t length, int srcpe ) const;
    virtual void put( void* dst, const void* src, size_t length, int srcpe ) const;

    virtual void writeAt( std::ostream& stream ) const;
    PartitionId mRank;
    PartitionId mSize;
    PGASAllocator * mAllocator;
private:
    static long mAllocated;
    static long mBarrierNum;
    long* pSync;
    void* pWork;
    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};
}

#endif // LAMA_GPIINTERFACE_HPP_
