/**
 * @file HArrayView.cpp
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
 * @brief Implementation and instantiation for template class HArrayView.
 * @author Lauretta Schubert
 * @date 22.08.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/HArrayView.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

namespace lama
{

/* HArrayView */

template<typename ValueType>
HArrayView<ValueType>::HArrayView( HArray<ValueType>& array )
                : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
    if( array.constFlag )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayView for a const HArray" )
    }
}

template<typename ValueType>
HArrayView<ValueType>::HArrayView( HArray<ValueType>& array, const IndexType offset, const IndexType size )
                : mArray( array ), mOffset( offset ), mSize( size )
{
    if( mOffset < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayView with a negative offset " << mOffset )
    }

    if( mSize < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayView with a negative size " << mSize )
    }

    if( mOffset + mSize > mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Could not create a HArrayView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }

    if( array.constFlag )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayView for a const HArray" )
    }
}

template<typename ValueType>
HArrayView<ValueType>::HArrayView( const HArrayView<ValueType>& other )
                : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename ValueType>
HArrayView<ValueType>::~HArrayView()
{
}

template<typename ValueType>
void HArrayView<ValueType>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename ValueType>
ValueType* HArrayView<ValueType>::get( const size_t index )
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
const ValueType* HArrayView<ValueType>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
int HArrayView<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename ValueType>
void HArrayView<ValueType>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename ValueType>
int HArrayView<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    return mArray.acquireWriteAccess( context, keepFlag );
}

template<typename ValueType>
int HArrayView<ValueType>::acquireWriteAccess()
{
    return mArray.acquireWriteAccess();
}

template<typename ValueType>
void HArrayView<ValueType>::releaseWriteAccess( const size_t index )
{
    mArray.releaseWriteAccess( index );
}

template<typename ValueType>
void HArrayView<ValueType>::clear( const size_t index )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Resizing a HArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.clear( index );
    mSize = 0;
}

template<typename ValueType>
void HArrayView<ValueType>::resize( const size_t index, const IndexType newSize )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Resizing a HArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.resize( index, newSize );
    mSize = mArray.size();
}

template<typename ValueType>
void HArrayView<ValueType>::reserve( const size_t index, const IndexType capacity, const bool copyFlag )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Calling reserve on a HArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.reserve( index, capacity, copyFlag );
}

template<typename ValueType>
IndexType HArrayView<ValueType>::capacity( const size_t index ) const
{
    IndexType capacityValue = mArray.capacity( index );

    if( mOffset != 0 || mSize != mArray.size() )
    {
        capacityValue = mSize; // no more available
    }

    return capacityValue;
}

/* HArrayConstView */

template<typename ValueType>
HArrayConstView<ValueType>::HArrayConstView( const HArrayConstView<ValueType>& other )
                : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename ValueType>
HArrayConstView<ValueType>::HArrayConstView( const HArrayView<ValueType>& view )
                : mArray( view.mArray ), mOffset( view.mOffset ), mSize( view.mSize )
{
}

template<typename ValueType>
HArrayConstView<ValueType>::HArrayConstView( const HArray<ValueType>& array )
                : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
}

template<typename ValueType>
HArrayConstView<ValueType>::HArrayConstView(
    const HArray<ValueType>& array,
    const IndexType offset,
    const IndexType size )
                : mArray( array ), mOffset( offset ), mSize( size )
{
    if( mOffset < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayConstView with a negative offset " << mOffset )
    }

    if( mSize < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a HArrayConstView with a negative size " << mSize )
    }

    if( mOffset + mSize > mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Could not create a HArrayConstView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }
}

template<typename ValueType>
HArrayConstView<ValueType>::~HArrayConstView()
{
}

template<typename ValueType>
void HArrayConstView<ValueType>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename ValueType>
const ValueType* HArrayConstView<ValueType>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
int HArrayConstView<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename ValueType>
void HArrayConstView<ValueType>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename ValueType>
bool HArrayConstView<ValueType>::operator==( const HArrayConstView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool HArrayConstView<ValueType>::operator!=( const HArrayConstView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool HArrayConstView<ValueType>::operator==( const HArrayView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool HArrayConstView<ValueType>::operator!=( const HArrayView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool HArrayView<ValueType>::operator==( const HArrayView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool HArrayView<ValueType>::operator!=( const HArrayView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool HArrayView<ValueType>::operator==( const HArrayConstView<ValueType>& other ) const
{
    return other == *this;
}

template<typename ValueType>
bool HArrayView<ValueType>::operator!=( const HArrayConstView<ValueType>& other ) const
{
    return !( *this == other );
}

// template instantiation for the supported data types

#define LAMA_ARRAY_VIEW_INSTANTIATE(z, I, _)                                   \
    template class COMMON_DLL_IMPORTEXPORT HArrayView< ARRAY_TYPE##I >;        \
    template class COMMON_DLL_IMPORTEXPORT HArrayConstView< ARRAY_TYPE##I >;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_ARRAY_VIEW_INSTANTIATE, _ )

#undef LAMA_ARRAY_VIEW_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
