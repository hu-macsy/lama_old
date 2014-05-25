/**
 * @file LAMAArrayView.cpp
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
 * @brief Implementation and instantiation for template class LAMAArrayView.
 * @author Lauretta Schubert
 * @date 22.08.2012
 * @since 1.0.0
 */

// hpp
#include <lama/LAMAArrayView.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace lama
{

/* LAMAArrayView */

template<typename T>
LAMAArrayView<T>::LAMAArrayView( LAMAArray<T>& array )
    : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
    if ( array.constFlag )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayView for a const LAMAArray" )
    }
}

template<typename T>
LAMAArrayView<T>::LAMAArrayView( LAMAArray<T>& array, const IndexType offset, const IndexType size )
    : mArray( array ), mOffset( offset ), mSize( size )
{
    if ( mOffset < 0 )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayView with a negative offset " << mOffset )
    }
    if ( mSize < 0 )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayView with a negative size " << mSize )
    }
    if ( mOffset + mSize > mArray.size() )
    {
        LAMA_THROWEXCEPTION(
            "Could not create a LAMAArrayView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }
    if ( array.constFlag )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayView for a const LAMAArray" )
    }
}

template<typename T>
LAMAArrayView<T>::LAMAArrayView( const LAMAArrayView<T>& other )
    : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename T>
LAMAArrayView<T>::~LAMAArrayView()
{
}

template<typename T>
void LAMAArrayView<T>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename T>
T* LAMAArrayView<T>::get( const size_t index )
{
    return mArray.get( index ) + mOffset;
}

template<typename T>
const T* LAMAArrayView<T>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename T>
int LAMAArrayView<T>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename T>
void LAMAArrayView<T>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename T>
int LAMAArrayView<T>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    return mArray.acquireWriteAccess( context, keepFlag );
}

template<typename T>
int LAMAArrayView<T>::acquireWriteAccess()
{
    return mArray.acquireWriteAccess();
}

template<typename T>
void LAMAArrayView<T>::releaseWriteAccess( const size_t index )
{
    mArray.releaseWriteAccess( index );
}

template<typename T>
void LAMAArrayView<T>::clear( const size_t index )
{
    if ( mOffset != 0 || mSize != mArray.size() )
    {
        LAMA_THROWEXCEPTION(
            "Resizing a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }
    mArray.clear( index );
    mSize = 0;
}

template<typename T>
void LAMAArrayView<T>::resize( const size_t index, const IndexType newSize )
{
    if ( mOffset != 0 || mSize != mArray.size() )
    {
        LAMA_THROWEXCEPTION(
            "Resizing a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }
    mArray.resize( index, newSize );
    mSize = mArray.size();
}

template<typename T>
void LAMAArrayView<T>::reserve( const size_t index, const IndexType capacity, const bool copyFlag )
{
    if ( mOffset != 0 || mSize != mArray.size() )
    {
        LAMA_THROWEXCEPTION(
            "Calling reserve on a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }
    mArray.reserve( index, capacity, copyFlag );
}

template<typename T>
IndexType LAMAArrayView<T>::capacity( const size_t index ) const
{
    IndexType capacityValue = mArray.capacity( index );

    if ( mOffset != 0 || mSize != mArray.size() )
    {
        capacityValue = mSize; // no more available
    }

    return capacityValue;
}

/* LAMAArrayConstView */

template<typename T>
LAMAArrayConstView<T>::LAMAArrayConstView( const LAMAArrayConstView<T>& other )
    : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename T>
LAMAArrayConstView<T>::LAMAArrayConstView( const LAMAArrayView<T>& view )
    : mArray( view.mArray ), mOffset( view.mOffset ), mSize( view.mSize )
{
}

template<typename T>
LAMAArrayConstView<T>::LAMAArrayConstView( const LAMAArray<T>& array )
    : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
}

template<typename T>
LAMAArrayConstView<T>::LAMAArrayConstView( const LAMAArray<T>& array, const IndexType offset, const IndexType size )
    : mArray( array ), mOffset( offset ), mSize( size )
{
    if ( mOffset < 0 )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayConstView with a negative offset " << mOffset )
    }
    if ( mSize < 0 )
    {
        LAMA_THROWEXCEPTION( "Could not create a LAMAArrayConstView with a negative size " << mSize )
    }
    if ( mOffset + mSize > mArray.size() )
    {
        LAMA_THROWEXCEPTION(
            "Could not create a LAMAArrayConstView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }
}

template<typename T>
LAMAArrayConstView<T>::~LAMAArrayConstView()
{
}

template<typename T>
void LAMAArrayConstView<T>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename T>
const T* LAMAArrayConstView<T>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename T>
int LAMAArrayConstView<T>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename T>
void LAMAArrayConstView<T>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename T>
bool LAMAArrayConstView<T>::operator==( const LAMAArrayConstView<T>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename T>
bool LAMAArrayConstView<T>::operator!=( const LAMAArrayConstView<T>& other ) const
{
    return !( *this == other );
}

template<typename T>
bool LAMAArrayConstView<T>::operator==( const LAMAArrayView<T>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename T>
bool LAMAArrayConstView<T>::operator!=( const LAMAArrayView<T>& other ) const
{
    return !( *this == other );
}

template<typename T>
bool LAMAArrayView<T>::operator==( const LAMAArrayView<T>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename T>
bool LAMAArrayView<T>::operator!=( const LAMAArrayView<T>& other ) const
{
    return !( *this == other );
}

template<typename T>
bool LAMAArrayView<T>::operator==( const LAMAArrayConstView<T>& other ) const
{
    return other == *this;
}

template<typename T>
bool LAMAArrayView<T>::operator!=( const LAMAArrayConstView<T>& other ) const
{
    return !( *this == other );
}

// template instantiation for the supported data types

#define LAMA_ARRAY_VIEW_INSTANTIATE(z, I, _)                                   \
   template class LAMA_DLL_IMPORTEXPORT LAMAArrayView< ARRAY_TYPE##I >;        \
   template class LAMA_DLL_IMPORTEXPORT LAMAArrayConstView< ARRAY_TYPE##I >;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_ARRAY_VIEW_INSTANTIATE, _ )

#undef LAMA_ARRAY_VIEW_INSTANTIATE

} /* namespace lama */
