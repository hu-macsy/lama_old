/**
 * @file PGASInterface.hpp
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
 * @brief PGASInterface.hpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */
#ifndef LAMA_PGASINTERFACE_HPP_
#define LAMA_PGASINTERFACE_HPP_

#include <lama/SyncToken.hpp>

#include <lama/pgas/PGASSyncToken.hpp>
#include <lama/pgas/PGASFunctor.hpp>

namespace lama
{

enum PGASCommunicationKind
{
    PGASget, PGASput
};

class PGASInterface: public Printable
{
//PAGSInterface
public:
    virtual void* allocate( size_t size ) const = 0;
    virtual void free( void* ptr, const size_t size ) const = 0;
    virtual void syncronizeAll() const = 0;
    virtual bool isPinned( const void* const ptr ) const = 0;
    virtual PGASSyncToken* getSyncToken( int arg1 ) const = 0;

    virtual PGASCommunicationKind getPreferredCommunicationKind() const = 0;

    virtual void get( void* dst, const void* src, size_t length, int srcpe ) const = 0;
    virtual void put( void* dst, const void* src, size_t length, int srcpe ) const = 0;

    virtual PartitionId getRank() const = 0;
    virtual PartitionId getSize() const = 0;

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

    virtual double max( const double val, const PartitionId root ) const;
    virtual float max( const float val, const PartitionId root ) const;
    virtual int max( const int val, const PartitionId root ) const;
    virtual size_t max( const size_t val, const PartitionId root ) const;

    virtual double min( const double val, const PartitionId root ) const;
    virtual float min( const float val, const PartitionId root ) const;
    virtual int min( const int val, const PartitionId root ) const;

    virtual double sum( const double val, const PartitionId root ) const;
    virtual float sum( const float val, const PartitionId root ) const;
    virtual int sum( const int val, const PartitionId root ) const;
    virtual size_t sum( const size_t val, const PartitionId root ) const;

    virtual double maxToAll( const double val ) const;
    virtual float maxToAll( const float val ) const;
    virtual int maxToAll( const int val ) const;
    virtual size_t maxToAll( const size_t val ) const;

    virtual double minToAll( const double val ) const;
    virtual float minToAll( const float val ) const;
    virtual int minToAll( const int val ) const;

    virtual double sumToAll( const double val ) const;
    virtual float sumToAll( const float val ) const;
    virtual int sumToAll( const int val ) const;
    virtual size_t sumToAll( const size_t val ) const;

    virtual void maxloc( double& d, int& loc, PartitionId root ) const;
    virtual void maxloc( float& d, int& loc, PartitionId root ) const;
    virtual void maxloc( int& d, int& loc, PartitionId root ) const;

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

    virtual void parallelReduction( PGASFunctor& reduction, PartitionId root ) const;
protected:
    virtual void writeAt( std::ostream& stream ) const = 0;
private:
    template<typename T>
    T maxToAllImpl( const T val ) const;

    template<typename T>
    T minToAllImpl( const T val ) const;

    template<typename T>
    T sumToAllImpl( const T val ) const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger );

//Singleton Pattern
public:
    static const PGASInterface* getInstance();
    PGASInterface();
    virtual ~PGASInterface();
private:
    static std::auto_ptr<PGASInterface> sInstance;
    static PGASInterface* init();
};

template<typename T>
T PGASInterface::maxToAllImpl( T val ) const
{
    syncronizeAll();
    T* res = static_cast<T*>( allocate( sizeof(T) ) );
    T* res1 = static_cast<T*>( allocate( sizeof(T) ) );
    syncronizeAll();
    *res = max( val, 0 );
    syncronizeAll();
    broadcast( res1, res, sizeof(T), 0 );
    T temp = *res1;
    free( res, sizeof(T) );
    free( res1, sizeof(T) );
    return temp;
}
template<typename T>
T PGASInterface::minToAllImpl( T val ) const
{
    T* res = static_cast<T*>( allocate( sizeof(T) ) );
    *res = min( val, 0 );
    syncronizeAll();
    broadcast( res, res, sizeof(T), 0 );
    T temp = *res;
    free( res, sizeof(T) );
    return temp;
}
template<typename T>
T PGASInterface::sumToAllImpl( T val ) const
{
    T* res = static_cast<T*>( allocate( sizeof(T) ) );
    *res = sum( val, 0 );
    syncronizeAll();
    broadcast( res, res, sizeof(T), 0 );
    T temp = *res;
    free( res, sizeof(T) );
    return temp;
}
}

#endif // LAMA_PGASINTERFACE_HPP_
