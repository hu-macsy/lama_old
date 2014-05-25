/**
 * @file CSRStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation and instantiation for template class ReadAccess.
 * @author Lauretta Schubert
 * @date 04.06.2011
 * @since 1.0.0
 */

// hpp
#include <lama/ReadAccess.hpp>

// boost
#include <boost/preprocessor.hpp>


namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ReadAccess<ValueType>::logger, "ReadAccess" )

/* --------------------------------------------------------------------------- */

template<typename T>
ReadAccess<T>::ReadAccess( const LAMAArray<ValueType>& array, ContextPtr context )
    : mArrayView( new LAMAArrayConstView<ValueType>(array,0,array.size() ) )
{
    LAMA_ASSERT_ERROR( context, "NULL context for read access");
    //LAMA_LOG_DEBUG(logger, "will acquire read access for " << *mArrayView << " at " << *context )
    mIndex = mArrayView->acquireReadAccess(context);
    //LAMA_LOG_TRACE(logger, "acquired read access for " << *mArrayView << " at " << *context );
}

template<typename T>
ReadAccess<T>::ReadAccess( const LAMAArrayView<ValueType>& view, ContextPtr context )
    : mArrayView( new LAMAArrayConstView<ValueType>( view ) )
{
    LAMA_ASSERT_ERROR( context, "NULL context for read access" )
    //LAMA_LOG_DEBUG(logger, "acquire read access for " << *mArrayView << " at " << *context);
    mIndex = mArrayView->acquireReadAccess( context );
    //LAMA_LOG_TRACE(logger, "acquired read access for " << *mArrayView << " at " << *context );
}

template<typename T>
ReadAccess<T>::ReadAccess( const LAMAArrayConstView<ValueType>& view, ContextPtr context )
    : mArrayView( new LAMAArrayConstView<ValueType>( view ) )
{
    LAMA_ASSERT_ERROR( context, "NULL context for read access" )
    //LAMA_LOG_DEBUG(logger, "acquire read access for " << *mArrayView << " at " << *context);
    mIndex = mArrayView->acquireReadAccess( context );
    //LAMA_LOG_TRACE(logger, "acquired read access for " << *mArrayView << " at " << *context );
}

template<typename T>
ReadAccess<T>::~ReadAccess()
{
    LAMA_LOG_TRACE( logger, "~ReadAccess: release" )
    release();
}

template<typename T>
const T* ReadAccess<T>::get() const
{
    if ( mArrayView == 0 )
    {
        LAMA_THROWEXCEPTION( "ReadAccess::get fails, has already been released." )
    }

    const T* ptr = mArrayView->get( mIndex );

    //LAMA_LOG_TRACE( logger, "mArray->get(" << mIndex<< ") with ptr = " << ptr )
    //*mArray->mContextData[mIndex].context << ", ptr = " << ptr )
    return ptr;
}

template<typename T>
void ReadAccess<T>::release()
{
    if ( mArrayView )
    {
        //LAMA_LOG_DEBUG(logger, "release read access for " << *mArrayView  );
        //LAMA_LOG_DEBUG(logger, "release read access for " << *mArray
        //          << " at " << *mArray->mContextData[mIndex].context);
        mArrayView->releaseReadAccess( mIndex );
        delete mArrayView;
    }
    else
    {
        //LAMA_LOG_TRACE(logger, "release read access for an already released array");
    }
    mArrayView = 0;
}

template<typename T>
void ReadAccess<T>::writeAt( std::ostream& stream ) const
{
    stream << "ReadAccess to ";
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

#define LAMA_READ_ACCESS_INSTANTIATE(z, I, _)                               \
template class LAMA_DLL_IMPORTEXPORT ReadAccess< ARRAY_TYPE##I >;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_READ_ACCESS_INSTANTIATE, _ )

#undef LAMA_READ_ACCESS_INSTANTIATE

} /* namespace lama */
