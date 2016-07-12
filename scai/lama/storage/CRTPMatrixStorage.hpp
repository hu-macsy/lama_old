/**
 * @file CRTPMatrixStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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
#include <scai/lama/mepr/CRTPMatrixStorageWrapper.hpp>

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
        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::setCSRDataImpl(
            static_cast<Derived*>( this ), numRows, numColumns,
            numValues, ia, ja, values, this->getContextPtr() );
    }

    /** Implementation of _MatrixStorage::setDIAData for derived class.
     *
     *  This routine requires that the derived class provides a corresponding
     *  routine setDIADataImpl that can deal with a typed version of values.
     */

    void setDIAData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::_HArray& values )
    {
        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::setDIADataImpl(
            static_cast<Derived*>( this ), numRows, numColumns, numValues, offsets, values, this->getContextPtr() );
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
        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::buildCSRDataImpl(
            static_cast<const Derived*>( this ), csrIA, csrJA, csrValues, this->getContextPtr() );
    }

    /** Get the i-th row of a storage as LAMA array. */

    void getRow( hmemo::_HArray& row, const IndexType irow ) const
    {
        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::getRowImpl( static_cast<const Derived*>( this ), row, irow );
    }

    void getDiagonal( hmemo::_HArray& diagonal ) const
    {
        if ( !this->hasDiagonalProperty() )
        {
            COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot get diagonal" )
        }

        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::getDiagonalImpl( static_cast<const Derived*>( this ), diagonal );
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

        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::setDiagonalVImpl( static_cast<Derived*>( this ), diagonal );
    }

    void scale( const ValueType value )
    {
        static_cast<Derived*>( this )->scaleImpl( value );
    }

    /** Polymorph implementation for MatrixStorage<ValueType>::scaleRows */

    void scaleRows( const hmemo::_HArray& values )
    {
        SCAI_ASSERT_EQUAL_ERROR( this->getNumRows(), values.size() )
        mepr::CRTPMatrixStorageWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::scaleRowsImpl( static_cast<Derived*>( this ), values );
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
