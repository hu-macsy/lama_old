/**
 * @file LAMAArrayView.cpp
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
 * @brief Implementation and instantiation for template class LAMAArrayView.
 * @author Lauretta Schubert
 * @date 22.08.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/LAMAArrayView.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace lama
{

/* LAMAArrayView */

template<typename ValueType>
LAMAArrayView<ValueType>::LAMAArrayView( LAMAArray<ValueType>& array )
                : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
    if( array.constFlag )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayView for a const LAMAArray" )
    }
}

template<typename ValueType>
LAMAArrayView<ValueType>::LAMAArrayView( LAMAArray<ValueType>& array, const IndexType offset, const IndexType size )
                : mArray( array ), mOffset( offset ), mSize( size )
{
    if( mOffset < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayView with a negative offset " << mOffset )
    }

    if( mSize < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayView with a negative size " << mSize )
    }

    if( mOffset + mSize > mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Could not create a LAMAArrayView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }

    if( array.constFlag )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayView for a const LAMAArray" )
    }
}

template<typename ValueType>
LAMAArrayView<ValueType>::LAMAArrayView( const LAMAArrayView<ValueType>& other )
                : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename ValueType>
LAMAArrayView<ValueType>::~LAMAArrayView()
{
}

template<typename ValueType>
void LAMAArrayView<ValueType>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename ValueType>
ValueType* LAMAArrayView<ValueType>::get( const size_t index )
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
const ValueType* LAMAArrayView<ValueType>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
int LAMAArrayView<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename ValueType>
void LAMAArrayView<ValueType>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename ValueType>
int LAMAArrayView<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    return mArray.acquireWriteAccess( context, keepFlag );
}

template<typename ValueType>
int LAMAArrayView<ValueType>::acquireWriteAccess()
{
    return mArray.acquireWriteAccess();
}

template<typename ValueType>
void LAMAArrayView<ValueType>::releaseWriteAccess( const size_t index )
{
    mArray.releaseWriteAccess( index );
}

template<typename ValueType>
void LAMAArrayView<ValueType>::clear( const size_t index )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Resizing a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.clear( index );
    mSize = 0;
}

template<typename ValueType>
void LAMAArrayView<ValueType>::resize( const size_t index, const IndexType newSize )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Resizing a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.resize( index, newSize );
    mSize = mArray.size();
}

template<typename ValueType>
void LAMAArrayView<ValueType>::reserve( const size_t index, const IndexType capacity, const bool copyFlag )
{
    if( mOffset != 0 || mSize != mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Calling reserve on a LAMAArrayView with 0 != offset = " << mOffset << " or a size not equal to the underlying array is not allowed. ( size = " << mSize << ", array size = " << mArray.size() );
    }

    mArray.reserve( index, capacity, copyFlag );
}

template<typename ValueType>
IndexType LAMAArrayView<ValueType>::capacity( const size_t index ) const
{
    IndexType capacityValue = mArray.capacity( index );

    if( mOffset != 0 || mSize != mArray.size() )
    {
        capacityValue = mSize; // no more available
    }

    return capacityValue;
}

/* LAMAArrayConstView */

template<typename ValueType>
LAMAArrayConstView<ValueType>::LAMAArrayConstView( const LAMAArrayConstView<ValueType>& other )
                : Printable(), mArray( other.mArray ), mOffset( other.mOffset ), mSize( other.mSize )
{
}

template<typename ValueType>
LAMAArrayConstView<ValueType>::LAMAArrayConstView( const LAMAArrayView<ValueType>& view )
                : mArray( view.mArray ), mOffset( view.mOffset ), mSize( view.mSize )
{
}

template<typename ValueType>
LAMAArrayConstView<ValueType>::LAMAArrayConstView( const LAMAArray<ValueType>& array )
                : mArray( array ), mOffset( 0 ), mSize( mArray.size() )
{
}

template<typename ValueType>
LAMAArrayConstView<ValueType>::LAMAArrayConstView(
    const LAMAArray<ValueType>& array,
    const IndexType offset,
    const IndexType size )
                : mArray( array ), mOffset( offset ), mSize( size )
{
    if( mOffset < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayConstView with a negative offset " << mOffset )
    }

    if( mSize < 0 )
    {
        COMMON_THROWEXCEPTION( "Could not create a LAMAArrayConstView with a negative size " << mSize )
    }

    if( mOffset + mSize > mArray.size() )
    {
        COMMON_THROWEXCEPTION(
                        "Could not create a LAMAArrayConstView with mOffset + mSize = " << mOffset + mSize << " because it exceeds the size of " << mArray )
    }
}

template<typename ValueType>
LAMAArrayConstView<ValueType>::~LAMAArrayConstView()
{
}

template<typename ValueType>
void LAMAArrayConstView<ValueType>::writeAt( std::ostream& stream ) const
{
    mArray.writeAt( stream );
}

template<typename ValueType>
const ValueType* LAMAArrayConstView<ValueType>::get( const size_t index ) const
{
    return mArray.get( index ) + mOffset;
}

template<typename ValueType>
int LAMAArrayConstView<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    return mArray.acquireReadAccess( context );
}

template<typename ValueType>
void LAMAArrayConstView<ValueType>::releaseReadAccess( const size_t index ) const
{
    mArray.releaseReadAccess( index );
}

template<typename ValueType>
bool LAMAArrayConstView<ValueType>::operator==( const LAMAArrayConstView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool LAMAArrayConstView<ValueType>::operator!=( const LAMAArrayConstView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool LAMAArrayConstView<ValueType>::operator==( const LAMAArrayView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool LAMAArrayConstView<ValueType>::operator!=( const LAMAArrayView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool LAMAArrayView<ValueType>::operator==( const LAMAArrayView<ValueType>& other ) const
{
    return &( other.mArray ) == &mArray;
}

template<typename ValueType>
bool LAMAArrayView<ValueType>::operator!=( const LAMAArrayView<ValueType>& other ) const
{
    return !( *this == other );
}

template<typename ValueType>
bool LAMAArrayView<ValueType>::operator==( const LAMAArrayConstView<ValueType>& other ) const
{
    return other == *this;
}

template<typename ValueType>
bool LAMAArrayView<ValueType>::operator!=( const LAMAArrayConstView<ValueType>& other ) const
{
    return !( *this == other );
}

// template instantiation for the supported data types

#define LAMA_ARRAY_VIEW_INSTANTIATE(z, I, _)                                   \
    template class COMMON_DLL_IMPORTEXPORT LAMAArrayView< ARRAY_TYPE##I >;        \
    template class COMMON_DLL_IMPORTEXPORT LAMAArrayConstView< ARRAY_TYPE##I >;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_ARRAY_VIEW_INSTANTIATE, _ )

#undef LAMA_ARRAY_VIEW_INSTANTIATE

} /* namespace lama */
