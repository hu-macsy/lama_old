/**
 * @file PGASCheckInterface.hpp
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
 * @brief PGASCheckInterface.hpp
 * @author Michael Drost
 * @date 07.02.2012
 * $Id$
 */
#ifndef LAMA_PGASCHECKINTERFACE_HPP_
#define LAMA_PGASCHECKINTERFACE_HPP_

#include <lama/pgas/PGASInterface.hpp>

namespace lama
{
class PGASCheckInterface: public PGASInterface
{
public:
    virtual void* allocate( size_t size ) const;
    virtual void free( void* ptr, const size_t size ) const;
    virtual void syncronizeAll() const;
    virtual bool isPinned( const void* const ptr ) const;
    virtual PGASSyncToken* getSyncToken( int arg1 ) const;

    virtual PGASCommunicationKind getPreferredCommunicationKind() const;

    virtual void get( void* dst, const void* src, size_t length, int srcpe ) const;
    virtual void put( void* dst, const void* src, size_t length, int srcpe ) const;

    virtual PartitionId getRank() const;
    virtual PartitionId getSize() const;

    virtual std::auto_ptr<SyncToken> getAsync( void* dst, const void* src, size_t length, int srcPE ) const;
    virtual std::auto_ptr<SyncToken> putAsync( void* dst, const void* src, size_t length, int srcPE ) const;

    virtual std::auto_ptr<SyncToken> shift(
        void* dst,
        const void* src,
        size_t size,
        PartitionId destRank,
        PartitionId srcRank ) const;

    virtual std::auto_ptr<SyncToken> broadcast( void* dst, const void* src, size_t length, int srcPE ) const;
    virtual std::auto_ptr<SyncToken> all2all( void* dst, const void* src, size_t elemSize ) const;

    virtual void swap( void* val, const size_t n, const PartitionId partner ) const;

    virtual void scatter( void* myvals, const size_t partSize, const PartitionId root, const void* allvals ) const;

    virtual void scatter(
        void* myvals,
        const size_t elemSize,
        const PartitionId root,
        const void* allvals,
        const IndexType sizes[] ) const;

    virtual void gather( void* allvals, const size_t partSize, const PartitionId root, const void* myvals ) const;

    virtual void gather(
        void* allvals,
        const size_t elemSize,
        const PartitionId root,
        const void* myvals,
        const IndexType sizes[] ) const;
    PGASCheckInterface( PGASInterface* interface );
    virtual ~PGASCheckInterface();
private:
    virtual void writeAt( std::ostream& stream ) const;
    std::auto_ptr<PGASInterface> mInterface;
};
}

#endif // LAMA_PGASCHECKINTERFACE_HPP_
