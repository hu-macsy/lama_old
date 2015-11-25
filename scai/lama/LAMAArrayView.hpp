/**
 * @file HArrayView.hpp
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
 * @brief HArrayView.hpp
 * @author Jiri Kraus
 * @date 05.08.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/lama/HArray.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class ReadAccess;

template<typename ValueType>
class WriteAccess;

template<typename ValueType>
class HArrayConstView;

/**
 * @brief HArrayView is a proxy which gives a view to a sub range of a HArray.
 *
 * @tparam ValueType is the type stored in this container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT HArrayView: public scai::common::Printable
{
    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;
    friend class HArrayConstView<ValueType> ;

public:

    /**
     * @brief Constructs a view into the whole passed HArray.
     *
     * @param[in,out] array the array to get a view into.
     */
    HArrayView( HArray<ValueType>& array );

    /**
     * @brief Constructs a view of the passed size into the passed HArray starting at the passed offset.
     *
     * @param[in,out]   array   the array to get a view into.
     * @param[in]       offset  the index where the view into array should start.
     * @param[in]       size    the lenght of the subrange of array.
     */
    HArrayView( HArray<ValueType>& array, const IndexType offset, const IndexType size );

    /**
     * @brief Takes a copy of the passed HArrayView.
     *
     * @param[in]   other   the HArrayView to take a copy from.
     */
    HArrayView( const HArrayView<ValueType>& other );

    /**
     * @brief destroys this HArrayView, the referenced HArray is not touched.
     */
    virtual ~HArrayView();

    /**
     * @brief Returns the size of this HArrayView.
     */
    inline IndexType size() const;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Checks if the HArray referenced by other is the same than the HArray reference by this.
     *
     * Checks if the the HArray referenced by other is the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the HArrayConstView to compare this with.
     * @return              if this and other are referencing the same HArray.
     */
    bool operator==( const HArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is not the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is not the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the HArrayConstView to compare this with.
     * @return              if this and other are referencing different HArrays.
     */
    bool operator!=( const HArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the HArrayView to compare this with.
     * @return              if this and other are referencing the same HArray.
     */
    bool operator==( const HArrayView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is not the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is not the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the HArrayView to compare this with.
     * @return              if this and other are referencing different HArrays.
     */
    bool operator!=( const HArrayView<ValueType>& other ) const;

private:
    HArrayView();
    HArrayView& operator=( const HArrayView<ValueType>& other );

    const ValueType* get( const size_t index ) const;

    ValueType* get( const size_t index );

    int acquireReadAccess( ContextPtr context ) const;

    void releaseReadAccess( const size_t index ) const;

    int acquireWriteAccess( ContextPtr context, bool keepFlag );

    int acquireWriteAccess();

    void releaseWriteAccess( const size_t index );

    void clear( const size_t index );

    void resize( const size_t index, const IndexType newSize );

    void reserve( const size_t index, const IndexType capacity, bool copyFlag );

    IndexType capacity( const size_t index ) const;

    HArray<ValueType>& mArray;
    const IndexType mOffset;
    IndexType mSize;
};

/**
 * @brief HArrayConstView is a proxy which gives a constant view to a sub range of a HArray.
 *
 * @tparam ValueType the value type for the elements of this.
 */
template<typename ValueType>
class HArrayConstView: public scai::common::Printable
{
    friend class ReadAccess<ValueType> ;

public:

    /**
     * @brief Takes a copy of the passed HArrayConstView.
     *
     * @param[in]   other   the HArrayConstView to take a copy from.
     */
    HArrayConstView( const HArrayConstView<ValueType>& other );

    /**
     * @brief Takes a copy of the passed HArrayView.
     *
     * @param[in]   view   the HArrayView to take a copy from.
     */
    HArrayConstView( const HArrayView<ValueType>& view );

    /**
     * @brief Constructs a view into the whole passed HArray.
     *
     * @param[in,out] array the array to get a view into.
     */
    HArrayConstView( const HArray<ValueType>& array );

    /**
     * @brief Constructs a view of the passed size into the passed HArray starting at the passed offset.
     *
     * @param[in,out]   array   the array to get a view into.
     * @param[in]       offset  the index where the view into array should start.
     * @param[in]       size    the lenght of the subrange of array.
     */
    HArrayConstView( const HArray<ValueType>& array, const IndexType offset, const IndexType size );

    /**
     * @brief destroys this HArrayView, the referenced HArray is not touched.
     */
    virtual ~HArrayConstView();

    /**
     * @brief Returns the size of this HArrayView.
     */
    inline IndexType size() const;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Checks if the HArray referenced by other is the same than the HArray reference by this.
     *
     * Checks if the the HArray referenced by other is the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the HArrayConstView to compare this with.
     * @return              if this and other are referencing the same HArray.
     */
    bool operator==( const HArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is not the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is not the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the HArrayConstView to compare this with.
     * @return              if this and other are referencing different HArrays.
     */
    bool operator!=( const HArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the HArrayView to compare this with.
     * @return              if this and other are referencing the same HArray.
     */
    bool operator==( const HArrayView<ValueType>& other ) const;

    /**
     * @brief Checks if the HArray referenced by other is not the same than the HArray reference by this.
     *
     * Checks if the HArray referenced by other is not the same than the HArray reference by this.
     * CAVEAT: It is only checked if the HArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the HArrayView to compare this with.
     * @return              if this and other are referencing different HArrays.
     */
    bool operator!=( const HArrayView<ValueType>& other ) const;

private:
    HArrayConstView();
    HArrayConstView& operator=( const HArrayConstView<ValueType>& other );

    const ValueType* get( const size_t index ) const;

    int acquireReadAccess( ContextPtr context ) const;

    void releaseReadAccess( const size_t index ) const;

    const HArray<ValueType>& mArray;
    const IndexType mOffset;
    const IndexType mSize;
};

template<typename ValueType>
inline IndexType HArrayView<ValueType>::size() const
{
    return mSize;
}

template<typename ValueType>
inline IndexType HArrayConstView<ValueType>::size() const
{
    return mSize;
}

} /* end namespace lama */

} /* end namespace scai */
