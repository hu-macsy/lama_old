/**
 * @file LAMAArrayView.hpp
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
 * @brief LAMAArrayView.hpp
 * @author Jiri Kraus
 * @date 05.08.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// others
#include <scai/lama/LAMAArray.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class ReadAccess;

template<typename ValueType>
class WriteAccess;

template<typename ValueType>
class LAMAArrayConstView;

/**
 * @brief LAMAArrayView is a proxy which gives a view to a sub range of a LAMAArray.
 *
 * @tparam ValueType is the type stored in this container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArrayView: public Printable
{
    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;
    friend class LAMAArrayConstView<ValueType> ;

public:

    /**
     * @brief Constructs a view into the whole passed LAMAArray.
     *
     * @param[in,out] array the array to get a view into.
     */
    LAMAArrayView( LAMAArray<ValueType>& array );

    /**
     * @brief Constructs a view of the passed size into the passed LAMAArray starting at the passed offset.
     *
     * @param[in,out]   array   the array to get a view into.
     * @param[in]       offset  the index where the view into array should start.
     * @param[in]       size    the lenght of the subrange of array.
     */
    LAMAArrayView( LAMAArray<ValueType>& array, const IndexType offset, const IndexType size );

    /**
     * @brief Takes a copy of the passed LAMAArrayView.
     *
     * @param[in]   other   the LAMAArrayView to take a copy from.
     */
    LAMAArrayView( const LAMAArrayView<ValueType>& other );

    /**
     * @brief destroys this LAMAArrayView, the referenced LAMAArray is not touched.
     */
    virtual ~LAMAArrayView();

    /**
     * @brief Returns the size of this LAMAArrayView.
     */
    inline IndexType size() const;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     *
     * Checks if the the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the LAMAArrayConstView to compare this with.
     * @return              if this and other are referencing the same LAMAArray.
     */
    bool operator==( const LAMAArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the LAMAArrayConstView to compare this with.
     * @return              if this and other are referencing different LAMAArrays.
     */
    bool operator!=( const LAMAArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the LAMAArrayView to compare this with.
     * @return              if this and other are referencing the same LAMAArray.
     */
    bool operator==( const LAMAArrayView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the LAMAArrayView to compare this with.
     * @return              if this and other are referencing different LAMAArrays.
     */
    bool operator!=( const LAMAArrayView<ValueType>& other ) const;

private:
    LAMAArrayView();
    LAMAArrayView& operator=( const LAMAArrayView<ValueType>& other );

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

    LAMAArray<ValueType>& mArray;
    const IndexType mOffset;
    IndexType mSize;
};

/**
 * @brief LAMAArrayConstView is a proxy which gives a constant view to a sub range of a LAMAArray.
 *
 * @tparam ValueType the value type for the elements of this.
 */
template<typename ValueType>
class LAMAArrayConstView: public Printable
{
    friend class ReadAccess<ValueType> ;

public:

    /**
     * @brief Takes a copy of the passed LAMAArrayConstView.
     *
     * @param[in]   other   the LAMAArrayConstView to take a copy from.
     */
    LAMAArrayConstView( const LAMAArrayConstView<ValueType>& other );

    /**
     * @brief Takes a copy of the passed LAMAArrayView.
     *
     * @param[in]   view   the LAMAArrayView to take a copy from.
     */
    LAMAArrayConstView( const LAMAArrayView<ValueType>& view );

    /**
     * @brief Constructs a view into the whole passed LAMAArray.
     *
     * @param[in,out] array the array to get a view into.
     */
    LAMAArrayConstView( const LAMAArray<ValueType>& array );

    /**
     * @brief Constructs a view of the passed size into the passed LAMAArray starting at the passed offset.
     *
     * @param[in,out]   array   the array to get a view into.
     * @param[in]       offset  the index where the view into array should start.
     * @param[in]       size    the lenght of the subrange of array.
     */
    LAMAArrayConstView( const LAMAArray<ValueType>& array, const IndexType offset, const IndexType size );

    /**
     * @brief destroys this LAMAArrayView, the referenced LAMAArray is not touched.
     */
    virtual ~LAMAArrayConstView();

    /**
     * @brief Returns the size of this LAMAArrayView.
     */
    inline IndexType size() const;

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     *
     * Checks if the the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the LAMAArrayConstView to compare this with.
     * @return              if this and other are referencing the same LAMAArray.
     */
    bool operator==( const LAMAArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the LAMAArrayConstView to compare this with.
     * @return              if this and other are referencing different LAMAArrays.
     */
    bool operator!=( const LAMAArrayConstView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are the same. Size and offset of the views are not checked for equality.
     *
     * @param[in]   other   the LAMAArrayView to compare this with.
     * @return              if this and other are referencing the same LAMAArray.
     */
    bool operator==( const LAMAArrayView<ValueType>& other ) const;

    /**
     * @brief Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     *
     * Checks if the LAMAArray referenced by other is not the same than the LAMAArray reference by this.
     * CAVEAT: It is only checked if the LAMAArrays are not the same. Size and offset of the views are not checked for inequality.
     *
     * @param[in]   other   the LAMAArrayView to compare this with.
     * @return              if this and other are referencing different LAMAArrays.
     */
    bool operator!=( const LAMAArrayView<ValueType>& other ) const;

private:
    LAMAArrayConstView();
    LAMAArrayConstView& operator=( const LAMAArrayConstView<ValueType>& other );

    const ValueType* get( const size_t index ) const;

    int acquireReadAccess( ContextPtr context ) const;

    void releaseReadAccess( const size_t index ) const;

    const LAMAArray<ValueType>& mArray;
    const IndexType mOffset;
    const IndexType mSize;
};

template<typename ValueType>
inline IndexType LAMAArrayView<ValueType>::size() const
{
    return mSize;
}

template<typename ValueType>
inline IndexType LAMAArrayConstView<ValueType>::size() const
{
    return mSize;
}

} /* end namespace lama */

} /* end namespace scai */
