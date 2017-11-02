/**
 * @file lama/mepr/CRTPMatrixStorageWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Wrapper for templated calls in CRTPMatrixStorage
 * @author Eric Schricker
 * @date 10.03.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

namespace scai
{

namespace lama
{

namespace mepr
{

/*
 * Forward declaration
 */
template<typename Derived, typename TList>
struct CRTPMatrixStorageWrapper;

/*
 * Termination
 */
template<typename Derived>
struct CRTPMatrixStorageWrapper<Derived, common::mepr::NullType>
{
    static void setCSRDataImpl(
        Derived*,
        const IndexType,
        const IndexType,
        const IndexType,
        const hmemo::HArray<IndexType>&,
        const hmemo::HArray<IndexType>&,
        const hmemo::_HArray&,
        hmemo::ContextPtr )
    {}

    static void setDIADataImpl(
        Derived*,
        const IndexType,
        const IndexType,
        const IndexType,
        const hmemo::HArray<IndexType>&,
        const hmemo::_HArray&,
        hmemo::ContextPtr )
    {}

    static void buildCSRDataImpl(
        const Derived*,
        hmemo::HArray<IndexType>&,
        hmemo::HArray<IndexType>&,
        hmemo::_HArray&,
        const hmemo::ContextPtr )
    {}

    static void getRowImpl(
        const Derived*,
        hmemo::_HArray&,
        const IndexType )
    {}

    static void setRowImpl(
        Derived*,
        const hmemo::_HArray&,
        const IndexType,
        const common::BinaryOp )
    {}

    static void setColumnImpl(
        Derived*,
        const hmemo::_HArray&,
        const IndexType,
        const common::BinaryOp )
    {}

    static void getColumnImpl(
        const Derived*,
        hmemo::_HArray&,
        const IndexType )
    {}

    static void getDiagonalImpl(
        const Derived*,
        hmemo::_HArray& )
    {}

    static void setDiagonalVImpl(
        Derived*,
        const hmemo::_HArray& )
    {}

    static void scaleRowsImpl(
        Derived*,
        const hmemo::_HArray& )
    {}
};

/*
 * Step n
 */
template<typename Derived, typename H, typename T>
struct CRTPMatrixStorageWrapper<Derived, common::mepr::TypeList<H, T> >
{
    static void setCSRDataImpl(
        Derived* obj,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values,
        hmemo::ContextPtr ctx )
    {
        if ( values.getValueType() == common::getScalarType<H>() )
        {
            obj->setCSRDataImpl( numRows, numColumns, numValues, ia, ja, reinterpret_cast<const hmemo::HArray<H>& >( values ), ctx );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setCSRDataImpl( obj, numRows, numColumns, numValues, ia, ja, values, ctx );
        }
    }

    static void setDIADataImpl(
        Derived* obj,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::_HArray& values,
        hmemo::ContextPtr ctx )
    {
        if ( values.getValueType() == common::getScalarType<H>() )
        {
            obj->setDIADataImpl( numRows, numColumns, numDiagonals, offsets, reinterpret_cast<const hmemo::HArray<H>& >( values ), ctx );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setDIADataImpl( obj, numRows, numColumns, numDiagonals, offsets, values, ctx );
        }
    }

    static void buildCSRDataImpl(
        const Derived* obj,
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::_HArray& csrValues,
        const hmemo::ContextPtr ctx )
    {
        if ( csrValues.getValueType() == common::getScalarType<H>() )
        {
            obj->buildCSR( csrIA, &csrJA, reinterpret_cast<hmemo::HArray<H>* >( &csrValues ), ctx );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::buildCSRDataImpl( obj, csrIA, csrJA, csrValues, ctx );
        }
    }

    static void getRowImpl(
        const Derived* obj,
        hmemo::_HArray& row,
        const IndexType irow )
    {
        if ( row.getValueType() == common::getScalarType<H>() )
        {
            obj->getRowImpl( reinterpret_cast<hmemo::HArray<H>& >( row ), irow );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::getRowImpl( obj, row, irow );
        }
    }

    static void setRowImpl(
        Derived* obj,
        const hmemo::_HArray& row,
        const IndexType i,
        const common::BinaryOp op )
    {
        if ( row.getValueType() == common::getScalarType<H>() )
        {
            obj->setRowImpl( reinterpret_cast<const hmemo::HArray<H>& >( row ), i, op );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setRowImpl( obj, row, i, op );
        }
    }

    static void getColumnImpl(
        const Derived* obj,
        hmemo::_HArray& column,
        const IndexType j )
    {
        if ( column.getValueType() == common::getScalarType<H>() )
        {
            obj->getColumnImpl( reinterpret_cast<hmemo::HArray<H>& >( column ), j );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::getColumnImpl( obj, column, j );
        }
    }

    static void setColumnImpl(
        Derived* obj,
        const hmemo::_HArray& column,
        const IndexType j,
        const common::BinaryOp op )
    {
        if ( column.getValueType() == common::getScalarType<H>() )
        {
            obj->setColumnImpl( reinterpret_cast<const hmemo::HArray<H>& >( column ), j, op );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setColumnImpl( obj, column, j, op );
        }
    }

    static void getDiagonalImpl(
        const Derived* obj,
        hmemo::_HArray& diagonal )
    {
        if ( diagonal.getValueType() == common::getScalarType<H>() )
        {
            obj->getDiagonalImpl( reinterpret_cast<hmemo::HArray<H>& >( diagonal ) );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::getDiagonalImpl( obj, diagonal );
        }
    }

    static void setDiagonalVImpl(
        Derived* obj,
        const hmemo::_HArray& diagonal )
    {
        if ( diagonal.getValueType() == common::getScalarType<H>() )
        {
            obj->setDiagonalImpl( reinterpret_cast<const hmemo::HArray<H>& >( diagonal ) );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setDiagonalVImpl( obj, diagonal );
        }
    }

    static void scaleRowsImpl(
        Derived* obj,
        const hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<H>() )
        {
            obj->scaleImpl( reinterpret_cast<const hmemo::HArray<H>& >( values ) );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::scaleRowsImpl( obj, values );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
