/**
 * @file CRTPMatrixStorage.hpp
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
 * @brief Extension of a derived matrix storage class with additional routines.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/mepr/CRTPWrapper.hpp>

namespace scai
{

namespace lama
{

/** This template class supports static polymorphism to define
 *  common routines for derived classes of MatrixStorage.
 *
 *  Therefore it uses the Curiously Recurring Template Pattern (CRTP)
 *  as a C++ idiom where a class X derived from this template is
 *  itself a template argument.
 *
 *  @tparam Derived is the derived class of MatrixStorage<ValueType>
 *  @tparam ValueType specifies the type of the matrix valuess
 *
 *  @todo: create, copy should also be defined here
 */

template<class Derived, typename ValueType>
class COMMON_DLL_IMPORTEXPORT CRTPMatrixStorage: public MatrixStorage<ValueType>
{
public:

    /** Constructor for a zero-matrix of a certain size. */

    CRTPMatrixStorage( const IndexType numRows, const IndexType numColumns )

        : MatrixStorage<ValueType>( numRows, numColumns )

    {
    }

    /** Constructor for an empty matrix. */

    CRTPMatrixStorage()

        : MatrixStorage<ValueType>( 0, 0 )

    {
    }

    MatrixStorageCreateKeyType getCreateValue() const
    {
        return Derived::createValue();
    }

    /** Implementation of _MatrixStorage::setCSRData for derived class.
     *
     *  This routine requires that the derived class provides a corresponding
     *  routine setCSRDataImpl that can deal with a typed version of values.
     */

    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values )
    {
        mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::setCSRDataImpl(
                        static_cast<Derived*>( this ), numRows, numColumns,
                        numValues, ia, ja, values, this->getContextPtr());
    }

    /** Implementation for _MatrixStorage::buildCSRSizes */

    void buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
    {
// The sizes will be available via buildCSR with NULL for ja, values

        hmemo::HArray<IndexType>* ja = NULL;
        hmemo::HArray<ValueType>* values = NULL;

        static_cast<const Derived*>( this )->buildCSR( ia, ja, values, this->getContextPtr() );
    }

    void buildCSRData( hmemo::HArray<IndexType>& csrIA, hmemo::HArray<IndexType>& csrJA, hmemo::_HArray& csrValues ) const
    {
         mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::buildCSRDataImpl(
                         static_cast<const Derived*>( this ), csrIA, csrJA, csrValues, this->getContextPtr() );
    }

    /** Get the i-th row of a storage as LAMA array. */

    void getRow( hmemo::_HArray& row, const IndexType irow ) const
    {
        mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::getRowImpl( static_cast<const Derived* const>( this ), row, irow );
    }

    void getDiagonal( hmemo::_HArray& diagonal ) const
    {
        if ( !this->hasDiagonalProperty() )
        {
            COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot get diagonal" )
        }

        mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::getDiagonalImpl( static_cast<const Derived* const>( this ), diagonal );
    }

    void setDiagonal( const ValueType value )
    {
        if ( !this->hasDiagonalProperty() )
        {
            COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" )
        }

        static_cast<Derived*>( this )->setDiagonalImpl( value );
    }

    void setDiagonalV( const hmemo::_HArray& diagonal )
    {
        IndexType numDiagonalElements = diagonal.size();

        if ( numDiagonalElements > this->getNumRows() || numDiagonalElements > this->getNumColumns() )
        {
            COMMON_THROWEXCEPTION( "Diagonal of size " << numDiagonalElements << " too large for matrix: " << *this )
        }

        if ( !this->hasDiagonalProperty() )
        {
            COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" )
        }

        mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::setDiagonalVImpl( static_cast<Derived*>( this ), diagonal );
    }

    void scale( const ValueType value )
    {
        static_cast<Derived*>( this )->scaleImpl( value );
    }

    /** Polymorph implementation for MatrixStorage<ValueType>::scaleRows */

    void scaleRows( const hmemo::_HArray& values )
    {
        SCAI_ASSERT_EQUAL_ERROR( this->getNumRows(), values.size() )

        mepr::CRTPWrapper<Derived, ARITHMETIC_HOST_LIST>::scaleRowsImpl( static_cast<Derived*>( this ), values );
    }

    /** Implementation of MatrixStorage::getTypeName for derived class. */

    virtual const char* getTypeName() const
    {
// each derived class provides static method to get the type name.

        return Derived::typeName();
    }
};

} /* end namespace lama */

} /* end namespace scai */
