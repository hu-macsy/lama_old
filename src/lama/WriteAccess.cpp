/**
 * @file CSRStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Implementation and instantiation for template class WriteAccess.
 * @author Lauretta Schubert
 * @date 22.08.2012
 * @since 1.0.0
 */

// hpp
#include <lama/WriteAccess.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, WriteAccess<ValueType>::logger, "WriteAccess" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( LAMAArray<ValueType>& view, ContextPtr context, const bool keep /* = true*/)
                : mArrayView( &view ), mIndex( std::numeric_limits<size_t>::max() )
{
    if( view.constFlag )
    {
        LAMA_THROWEXCEPTION( "write on const array not allowed" )
    }

    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArrayView << " at " << *context << ", keep = " << keep )
    mIndex = mArrayView->acquireWriteAccess( context, keep );
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess(
    LAMAArray<ValueType>& array,
    ContextPtr context,
    const IndexType size,
    const bool keep /* = false */)

                : mArrayView( &array ), mIndex( std::numeric_limits<size_t>::max() )

{
    if( array.constFlag )
    {
        LAMA_THROWEXCEPTION( "write on const array not allowed" )
    }

    LAMA_LOG_DEBUG( logger,
                    "acquire write access for " << *mArrayView << " at " << *context << ", size = " << size << ", keep = " << keep )
    mIndex = mArrayView->acquireWriteAccess( context, keep );

    if( !keep )
    {
        mArrayView->clear( mIndex );
    }

    mArrayView->resize( mIndex, size );
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( LAMAArray<ValueType>& array )
                : mArrayView( &array ), mIndex( std::numeric_limits<size_t>::max() )
{
    if( array.constFlag )
    {
        LAMA_THROWEXCEPTION( "write on const array not allowed" )
    }

    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArrayView << " at first valid context " )
    mIndex = mArrayView->acquireWriteAccess();
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename ValueType>
WriteAccess<ValueType>::~WriteAccess()
{
    LAMA_LOG_TRACE( logger, "~WriteAccess: release" )
    release();
}

template<typename ValueType>
ValueType* WriteAccess<ValueType>::get()
{
    if( !mArrayView )
    {
        LAMA_THROWEXCEPTION( "illegal get(): access has already been released." )
    }

    LAMA_LOG_TRACE( logger, "mData = " << mData )
    return mData;
}

template<typename ValueType>
void WriteAccess<ValueType>::clear()
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )
    mArrayView->clear( mIndex );
    mData = 0;
    LAMA_LOG_DEBUG( logger, "cleared " << *mArrayView )
}

template<typename ValueType>
void WriteAccess<ValueType>::resize( const IndexType newSize )
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )
    // do not log before check of mArrayView
    LAMA_LOG_DEBUG( logger, "resize " << *mArrayView << " to new size " << newSize )
    mArrayView->resize( mIndex, newSize );
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData )
}

template<typename ValueType>
void WriteAccess<ValueType>::reserve( const IndexType capacity )
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )
    LAMA_LOG_DEBUG( logger, "reserve " << *mArrayView << " to new capacity " << capacity )
    mArrayView->reserve( mIndex, capacity, true ); // copy = true for old data
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData )
}

template<typename ValueType>
IndexType WriteAccess<ValueType>::capacity() const
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )
    return mArrayView->capacity( mIndex );
}

template<typename ValueType>
void WriteAccess<ValueType>::release()
{
    if( mArrayView )
    {
        LAMA_LOG_DEBUG( logger, "release write access for " << *mArrayView )
        //LOG_DEBUG(logger, "release write access for " << *mArray << " at " << *mArray->mContextData[mIndex].context);
        mArrayView->releaseWriteAccess( mIndex );
        //delete mArrayView;
    }

    else
    {
        LAMA_LOG_DEBUG( logger, "release write access for an already released LAMAArray" )
    }

    mArrayView = 0;
    mData = 0;
    mIndex = std::numeric_limits<size_t>::max();
}

template<typename ValueType>
void WriteAccess<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "WriteAccess to ";

    if( mArrayView )
    {
        stream << *mArrayView;
    }

    else
    {
        stream << "already releases array view.";
    }
}

// template instantiation for the supported data types

#define LAMA_WRITE_ACCESS_INSTANTIATE(z, I, _)                           \
    template class COMMON_DLL_IMPORTEXPORT WriteAccess< ARRAY_TYPE##I >;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_WRITE_ACCESS_INSTANTIATE, _ )

#undef LAMA_WRITE_ACCESS_INSTANTIATE

} /* namespace lama */
