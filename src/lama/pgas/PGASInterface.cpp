/**
 * @file PGASInterface.cpp
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
 * @brief PGASInterface.cpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASInterface.hpp>
#include <lama/pgas/PGASNoInterface.hpp>
#include <lama/pgas/PGASCheckInterface.hpp>
#include <lama/pgas/functor/PGASSumFunctor.hpp>
#include <lama/pgas/functor/PGASMinFunctor.hpp>
#include <lama/pgas/functor/PGASMaxFunctor.hpp>
#include <lama/pgas/functor/PGASMaxLocFunctor.hpp>

#include <lama/NoSyncToken.hpp>

#ifdef LAMA_BUILD_GPI
#include <lama/pgas/GPIInterface.hpp>
#else
#ifdef LAMA_BUILD_OSHMEM
#include <lama/pgas/OpenShMemInterface.hpp>
#endif
#endif
namespace lama
{
LAMA_LOG_DEF_LOGGER( PGASInterface::logger, "PGASInterface" );

std::auto_ptr<SyncToken> PGASInterface::all2all( void *dst, const void *src, size_t elemSize ) const
{
    std::auto_ptr<SyncToken> token( new PGASSyncToken() );
    //Split this up so not all Nodes access the same Node at a Time
    const void* srcptr = static_cast<const char*>( src ) + ( getRank() * elemSize );

    //Copy the local Part
    void* dstptr = static_cast<char*>( dst ) + ( getRank() * elemSize );
    memcpy( dstptr, srcptr, sizeof( elemSize ) );

    for( int i = getRank() + 1; i < getSize(); ++i )
    {
        void* dstptr = static_cast<char*>( dst ) + i * elemSize;
        token->pushSyncToken( getAsync( dstptr, srcptr, elemSize, i ) );
    }
    for( int i = 0; i < getRank(); ++i )
    {
        void* dstptr = static_cast<char*>( dst ) + i * elemSize;
        token->pushSyncToken( getAsync( dstptr, srcptr, elemSize, i ) );
    }
    return token;
}

std::auto_ptr<SyncToken> PGASInterface::shift(
    void *dst,
    const void *src,
    size_t size,
    PartitionId destRank,
    PartitionId srcRank ) const
{
    LAMA_ASSERT_ERROR( srcRank != getRank(), "source must not be this partition" );
    LAMA_ASSERT_ERROR( destRank != getRank(), "dest must not be this partition" );
    // need an PGAS communicator with 2 requests, no clean up needed
    std::auto_ptr<SyncToken> pSyncToken( new PGASSyncToken() );
    if( getPreferredCommunicationKind() == PGASget )
    {
        pSyncToken->pushSyncToken( getAsync( dst, src, size, srcRank ) );
    }
    else
    {
        pSyncToken->pushSyncToken( putAsync( dst, src, size, destRank ) );
    }
    return pSyncToken;
}

size_t PGASInterface::sum( const size_t val, const PartitionId root ) const
{
    PGASSumFunctor<size_t> sumFunctor( val );
    parallelReduction( sumFunctor, root );
    return sumFunctor.getResult();
}

double PGASInterface::sum( const double val, const PartitionId root ) const
{
    PGASSumFunctor<double> sumFunctor( val );
    parallelReduction( sumFunctor, root );
    return sumFunctor.getResult();
}

float PGASInterface::sum( const float val, const PartitionId root ) const
{
    PGASSumFunctor<float> sumFunctor( val );
    parallelReduction( sumFunctor, root );
    return sumFunctor.getResult();
}

int PGASInterface::sum( const int val, const PartitionId root ) const
{
    PGASSumFunctor<int> sumFunctor( val );
    parallelReduction( sumFunctor, root );
    return sumFunctor.getResult();
}

double PGASInterface::max( const double val, const PartitionId root ) const
{
    PGASMaxFunctor<double> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

float PGASInterface::max( const float val, const PartitionId root ) const
{
    PGASMaxFunctor<float> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

int PGASInterface::max( const int val, const PartitionId root ) const
{
    PGASMaxFunctor<int> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}
size_t PGASInterface::max( const size_t val, const PartitionId root ) const
{
    PGASMaxFunctor<size_t> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

double PGASInterface::min( const double val, const PartitionId root ) const
{
    PGASMinFunctor<double> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

float PGASInterface::min( const float val, const PartitionId root ) const
{
    PGASMinFunctor<float> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

int PGASInterface::min( const int val, const PartitionId root ) const
{
    PGASMinFunctor<int> functor( val );
    parallelReduction( functor, root );
    return functor.getResult();
}

void PGASInterface::parallelReduction( PGASFunctor & reduction, PartitionId root ) const
{
    int vRank = ( getRank() + getSize() - root ) % getSize();
    syncronizeAll();
    if( getSize() == 1 )
    {
        return;
    }

    int size = getSize();
    for( int i = 0; i < getSize(); ++i ) //Emergency Clause so we never get a infinite loop
    {
        int carry = size % 2;
        size /= 2;
        reduction.iteration( ( getRank() + size + carry ) % getSize(), vRank < size );
        size += carry;
        if( size <= 1 )
        {
            break;
        }
    }
    syncronizeAll();
}

double PGASInterface::maxToAll( const double val ) const
{
    return maxToAllImpl( val );
}

float PGASInterface::maxToAll( const float val ) const
{
    return maxToAllImpl( val );
}

int PGASInterface::maxToAll( const int val ) const
{
    return maxToAllImpl( val );
}

size_t PGASInterface::maxToAll( const size_t val ) const
{
    return maxToAllImpl( val );
}

double PGASInterface::minToAll( const double val ) const
{
    return minToAllImpl( val );
}

float PGASInterface::minToAll( const float val ) const
{
    return minToAllImpl( val );
}

int PGASInterface::minToAll( const int val ) const
{
    return minToAllImpl( val );
}

double PGASInterface::sumToAll( const double val ) const
{
    return sumToAllImpl( val );
}

float PGASInterface::sumToAll( const float val ) const
{
    return sumToAllImpl( val );
}

int PGASInterface::sumToAll( const int val ) const
{
    return sumToAllImpl( val );
}

size_t PGASInterface::sumToAll( const size_t val ) const
{
    return sumToAllImpl( val );
}

void PGASInterface::maxloc( double & d, int & loc, PartitionId root ) const
{
    PGASMaxLocFunctor<double> functor( d );
    parallelReduction( functor, root );
    if( getRank() == root )
    {
        d = functor.getResult();
        loc = functor.getLoc();
    }
}

void PGASInterface::maxloc( float & d, int & loc, PartitionId root ) const
{
    PGASMaxLocFunctor<float> functor( d );
    parallelReduction( functor, root );
    if( getRank() == root )
    {
        d = functor.getResult();
        loc = functor.getLoc();
    }
}

void PGASInterface::maxloc( int & d, int & loc, PartitionId root ) const
{
    PGASMaxLocFunctor<int> functor( d );
    parallelReduction( functor, root );
    if( getRank() == root )
    {
        d = functor.getResult();
        loc = functor.getLoc();
    }
}

std::auto_ptr<SyncToken> PGASInterface::getAsync( void *dst, const void *src, size_t length, int srcPE ) const
{
    get( dst, src, length, srcPE );
    return std::auto_ptr<SyncToken>( new NoSyncToken() );
}

std::auto_ptr<SyncToken> PGASInterface::putAsync( void *dst, const void *src, size_t length, int srcPE ) const
{
    put( dst, src, length, srcPE );
    return std::auto_ptr<SyncToken>( new NoSyncToken() );
}

std::auto_ptr<SyncToken> PGASInterface::broadcast( void *dst, const void *src, size_t length, int srcPE ) const
{
    if( getPreferredCommunicationKind() == PGASget )
    {
        if( getRank() != srcPE )
        {
            get( dst, src, length, srcPE );
        }
    }
    else
    {
        for( int i = 0; i < getSize(); i++ )
        {
            if( getRank() == i )
            {
                continue;
            }
            put( dst, src, length, i );
        }

    }

    if( getRank() == srcPE && src != dst )
    {
        memcpy( dst, src, length );
    }
    return std::auto_ptr<SyncToken>( new NoSyncToken() );
}

void PGASInterface::swap( void *val, const size_t n, const PartitionId partner ) const
{
    void *temp = allocate( n );
    if( getPreferredCommunicationKind() == PGASget )
    {
        memcpy( temp, val, n );
        syncronizeAll();
        get( val, temp, n, partner );
    }
    else
    {
        put( temp, val, n, partner );
        syncronizeAll();
        memcpy( val, temp, n );

    }
    free( temp, n );
}

void PGASInterface::scatter( void *myvals, const size_t partSize, const PartitionId root, const void *allvals ) const
{
    if( getPreferredCommunicationKind() == PGASget )
    {
        if( root != getRank() ) {
            get( myvals, static_cast<const char*>( allvals ) + getRank() * partSize, partSize, root );
        }

    }
    else
    {
        for( int i = 0; i < getSize(); i++ )
        {
            if( getRank() == i )
            {
                continue;
            }
            if( getRank() == root ) {
                put( myvals, static_cast<const char*>( allvals ) + i * partSize, partSize, i );
            }

        }
    }

    if( root == getRank() ) {
        memcpy( myvals, static_cast<const char*>( allvals ) + getRank() * partSize, partSize );
    }
}

void PGASInterface::scatter(
    void *myvals,
    const size_t elemSize,
    const PartitionId root,
    const void *allvals,
    const IndexType sizes[] ) const
{
    if( getPreferredCommunicationKind() == PGASget )
    {
        int offset = 0;
        for( int i = 0; i < getRank(); ++i )
        {
            offset += sizes[i];
        }
        char* srcptr = static_cast<char*>( const_cast<void*>( allvals ) );
        srcptr += offset * elemSize;
        get( myvals, srcptr, sizes[getRank()] * elemSize, root );
    }
    else
    {
        int offset = 0;
        for( int i = 0; i < getSize(); ++i )
        {
            offset += sizes[i];
            char* srcptr = static_cast<char*>( const_cast<void*>( allvals ) );
            srcptr += offset * elemSize;
            if( getRank() == root )
            {
                put( myvals, srcptr, sizes[i] * elemSize, i );
            }

        }
    }
}

void PGASInterface::gather( void *allvals, const size_t partSize, const PartitionId root, const void *myvals ) const
{
    if( getPreferredCommunicationKind() == PGASput )
    {
        if( getRank() != root ) {
            put( static_cast<char*>( allvals ) + ( getRank() * partSize ), myvals, partSize, root );
        }

    }
    else
    {
        if( getRank() == root )
        {
            for( int i = 0; i < getSize(); i++ )
            {
                if( getRank() == i )
                {
                    continue;
                }
                get( static_cast<char*>( allvals ) + ( i * partSize ), myvals, partSize, i );
            }
        }

    }

    if( getRank() == root ) {
        memcpy( static_cast<char*>( allvals ) + getRank() * partSize, myvals, partSize );
    }
}

void PGASInterface::gather(
    void *allvals,
    const size_t elemSize,
    const PartitionId root,
    const void *myvals,
    const IndexType sizes[] ) const
{
    if( getPreferredCommunicationKind() == PGASput )
    {
        int offset = 0;
        for( int i = 0; i < getRank(); ++i )
        {
            offset += sizes[i];
        }
        if( getRank() == root )
        {
            memcpy( static_cast<char*>( allvals ) + offset * elemSize, myvals, sizes[getRank()] * elemSize );
        }
        else
        {
            put( static_cast<char*>( allvals ) + offset * elemSize, myvals, sizes[getRank()] * elemSize, root );
        }
    }
    else
    {
        int offset = 0;
        if( getRank() == root )
        {
            for( int i = 0; i < getSize(); ++i )
            {
                if( getRank() == i )
                {
                    memcpy( static_cast<char*>( allvals ) + offset * elemSize, myvals, sizes[getRank()] * elemSize );
                    continue;
                }
                get( static_cast<char*>( allvals ) + offset * elemSize, myvals, sizes[i] * elemSize, i );
                offset += sizes[i];
            }
        }
    }

}
PGASInterface::PGASInterface()
{
    LAMA_LOG_DEBUG( logger, "PGASInterface()" );
}

PGASInterface::~PGASInterface()
{
    LAMA_LOG_DEBUG( logger, "~PGASInterface()" );
}

PGASInterface *PGASInterface::init()
{
#ifdef LAMA_BUILD_GPI
    PGASInterface *temp = new GPIInterface();
#else
#ifdef LAMA_BUILD_OSHMEM
    PGASInterface *temp = new OpenShMemInterface();
#else
    PGASInterface *temp = new PGASNoInterface();
#endif
#endif
    LAMA_LOG_INFO( logger, "PGASInterface initialized: " << *temp );
#ifdef LAMA_ASSERT_LEVEL_DEBUG
    return new PGASCheckInterface(temp);
#else
    return temp;
#endif
}

const PGASInterface* PGASInterface::getInstance()
{
    if( dynamic_cast<PGASNoInterface*>( sInstance.get() ) != NULL )
    {
        sInstance.reset( init() );
    }
    return sInstance.get();
}

std::auto_ptr<PGASInterface> PGASInterface::sInstance( new PGASNoInterface() );

}

