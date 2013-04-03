/**
 * @file CSRStorage.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @author lschubert
 * @date 22.08.2012
 * $Id$
 */

// hpp
#include <lama/WriteAccess.hpp>

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, WriteAccess<ValueType>::logger, "WriteAccess" )

/* --------------------------------------------------------------------------- */

template<typename T>
WriteAccess<T>::WriteAccess( LAMAArray<ValueType>& array,
                             ContextPtr context )
    : mArrayView( new LAMAArrayView<ValueType>( array, 0, array.size() ) ), mIndex(std::numeric_limits<size_t>::max())
{
    const bool keep = true;

    LAMA_LOG_DEBUG(logger, "acquire write access for " << *mArrayView << " at "
                   << *context << ", keep = " << keep);

    mIndex = mArrayView->acquireWriteAccess( context, keep );

    mData = mArrayView->get( mIndex );

    LAMA_LOG_TRACE(logger, "mData = " << mData << ", mIndex = " << mIndex);
}

template<typename T>
WriteAccess<T>::WriteAccess(
    LAMAArray<ValueType>& array,
    ContextPtr context,
    const IndexType size,
    const bool keep /* = false */)

    : mArrayView( new LAMAArrayView<ValueType>( array, 0, array.size() ) ), mIndex(
        std::numeric_limits<size_t>::max() )
{
    LAMA_LOG_DEBUG( logger,
                    "acquire write access for " << *mArrayView << " at " << *context << ", size = " << size << ", keep = " << keep )

    mIndex = mArrayView->acquireWriteAccess( context, keep );

    if ( !keep )
    {
        mArrayView->clear( mIndex );
    }

    mArrayView->resize( mIndex, size );

    mData = mArrayView->get( mIndex );

    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename T>
WriteAccess<T>::WriteAccess( LAMAArrayView<ValueType>& view, ContextPtr context, const bool keep /* = true*/)
    : mArrayView( new LAMAArrayView<ValueType>( view ) ), mIndex( std::numeric_limits<size_t>::max() )
{
    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArrayView << " at " << *context << ", keep = " << keep )

    mIndex = mArrayView->acquireWriteAccess( context, keep );

    mData = mArrayView->get( mIndex );

    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename T>
WriteAccess<T>::WriteAccess(
    LAMAArrayView<ValueType>& view,
    ContextPtr context,
    const IndexType size,
    const bool keep /* = false */)

    : mArrayView( new LAMAArrayView<ValueType>( view ) ), mIndex( std::numeric_limits<size_t>::max() )
{
    LAMA_LOG_DEBUG( logger,
                    "acquire write access for " << *mArrayView << " at " << *context << ", size = " << size << ", keep = " << keep )

    mIndex = mArrayView->acquireWriteAccess( context, keep );

    if ( !keep )
    {
        mArrayView->clear( mIndex );
    }

    mArrayView->resize( mIndex, size );

    mData = mArrayView->get( mIndex );

    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename T>
WriteAccess<T>::WriteAccess( LAMAArray<ValueType>& array )
    : mArrayView( new LAMAArrayView<ValueType>( array, 0, array.size() ) ), mIndex(
        std::numeric_limits<size_t>::max() )
{
    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArrayView << " at first valid context " )

    mIndex = mArrayView->acquireWriteAccess();

    mData = mArrayView->get( mIndex );

    LAMA_LOG_TRACE( logger, "mData = " << mData << ", mIndex = " << mIndex )
}

template<typename T>
WriteAccess<T>::~WriteAccess()
{
    LAMA_LOG_TRACE( logger, "~WriteAccess: release" )
    release();
}

template<typename T>
T* WriteAccess<T>::get()
{
    if ( !mArrayView )
    {
        LAMA_THROWEXCEPTION( "illegal get(): access has already been released." )
    }

    LAMA_LOG_TRACE( logger, "mData = " << mData )

    return mData;
}

template<typename T>
void WriteAccess<T>::clear()
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )

    mArrayView->clear( mIndex );
    mData = 0;
    LAMA_LOG_DEBUG( logger, "cleared " << *mArrayView )
}

template<typename T>
void WriteAccess<T>::resize( const IndexType newSize )
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )

    // do not log before check of mArrayView

    LAMA_LOG_DEBUG( logger, "resize " << *mArrayView << " to new size " << newSize )

    mArrayView->resize( mIndex, newSize );
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData )
}

template<typename T>
void WriteAccess<T>::reserve( const IndexType capacity )
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )

    LAMA_LOG_DEBUG( logger, "reserve " << *mArrayView << " to new capacity " << capacity )

    mArrayView->reserve( mIndex, capacity, true ); // copy = true for old data
    mData = mArrayView->get( mIndex );
    LAMA_LOG_TRACE( logger, "mData = " << mData )
}

template<typename T>
IndexType WriteAccess<T>::capacity() const
{
    LAMA_ASSERT_ERROR( mArrayView, "WriteAccess has already been released." )

    return mArrayView->capacity( mIndex );
}

template<typename T>
void WriteAccess<T>::release()
{
    if ( mArrayView )
    {
        LAMA_LOG_DEBUG( logger, "release write access for " << *mArrayView )
        //LOG_DEBUG(logger, "release write access for " << *mArray << " at " << *mArray->mContextData[mIndex].context);
        mArrayView->releaseWriteAccess( mIndex );
        delete mArrayView;
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "release write access for an already released LAMAArrayView" )
    }
    mArrayView = 0;
    mData = 0;
    mIndex = std::numeric_limits<size_t>::max();
}

template<typename T>
void WriteAccess<T>::writeAt( std::ostream& stream ) const
{
    stream << "WriteAccess to ";
    if ( mArrayView )
    {
        stream << *mArrayView;
    }
    else
    {
        stream << "already releases array view.";
    }
}

// template instantiation for the supported data types

template class LAMA_DLL_IMPORTEXPORT WriteAccess<IndexType> ;
template class LAMA_DLL_IMPORTEXPORT WriteAccess<float> ;
template class LAMA_DLL_IMPORTEXPORT WriteAccess<double> ;

} /* namespace lama */
