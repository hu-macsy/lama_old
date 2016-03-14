/*
 * CRTPStorage.hpp
 *
 *  Created on: Mar 10, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace lama {

namespace mepr {

template<typename Derived, typename TList>
struct CRTPMatrixStorageWrapper;

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

    static void buildCSRDataImpl(
        const Derived * const,
        hmemo::HArray<IndexType>&,
        hmemo::HArray<IndexType>&,
        hmemo::_HArray&,
        const hmemo::ContextPtr )
    {}

    static void getRowImpl(
        const Derived * const,
        hmemo::_HArray&,
        const IndexType )
    { }

    static void getDiagonalImpl(
        const Derived* const,
        hmemo::_HArray& )
    { }

    static void setDiagonalVImpl(
        Derived*,
        const hmemo::_HArray& )
    { }

    static void scaleRowsImpl(
        Derived*,
        const hmemo::_HArray& )
    { }
};

template<typename Derived, typename H, typename T>
struct CRTPMatrixStorageWrapper<Derived, common::mepr::TypeList<H,T> >
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
        if( values.getValueType() == common::getScalarType<H>() )
        {
            obj->setCSRDataImpl( numRows, numColumns, numValues, ia, ja, reinterpret_cast<const hmemo::HArray<H>& >( values ), ctx );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::setCSRDataImpl( obj, numRows, numColumns, numValues, ia, ja, values, ctx );
        }
    }

    static void buildCSRDataImpl(
            const Derived* const obj,
            hmemo::HArray<IndexType>& csrIA,
            hmemo::HArray<IndexType>& csrJA,
            hmemo::_HArray& csrValues,
            const hmemo::ContextPtr ctx )
    {
        if( csrValues.getValueType() == common::getScalarType<H>() )
        {
            obj->buildCSR( csrIA, &csrJA, reinterpret_cast<hmemo::HArray<H>* >( &csrValues ), ctx );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::buildCSRDataImpl( obj, csrIA, csrJA, csrValues, ctx );
        }
    }

    static void getRowImpl(
            const Derived * const obj,
            hmemo::_HArray& row,
            const IndexType irow )
    {
        if( row.getValueType() == common::getScalarType<H>() )
        {
            obj->getRowImpl( reinterpret_cast<hmemo::HArray<H>& >( row ), irow );
        }
        else
        {
            CRTPMatrixStorageWrapper<Derived, T>::getRowImpl( obj, row, irow );
        }
    }

    static void getDiagonalImpl(
        const Derived* const obj,
        hmemo::_HArray& diagonal )
    {
        if( diagonal.getValueType() == common::getScalarType<H>() )
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
        if( diagonal.getValueType() == common::getScalarType<H>() )
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
        if( values.getValueType() == common::getScalarType<H>() )
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
