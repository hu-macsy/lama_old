/**
 * @file lama/mepr/StorageWrapper.hpp
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
 * @brief Wrapper for templated calls in derived Storage classes
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
struct StorageWrapper;

/*
 * Termination
 */
template<typename Derived>
struct StorageWrapper<Derived, common::mepr::NullType>
{
    static void setCSRDataImpl(
        Derived*,
        const IndexType,
        const IndexType,
        const hmemo::HArray<IndexType>&,
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

    static void assignImpl(
        Derived* obj,
        const _MatrixStorage& otherStorage )
    {
        COMMON_THROWEXCEPTION( "assign to " << obj << ": unsupported type = " << otherStorage.getValueType() );
    }
};

/*
 * Step n
 */
template<typename Derived, typename H, typename T>
struct StorageWrapper<Derived, common::mepr::TypeList<H, T> >
{
    static void setCSRDataImpl(
        Derived* obj,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values,
        hmemo::ContextPtr ctx )
    {
        if ( values.getValueType() == common::getScalarType<H>() )
        {
            obj->setCSRDataImpl( numRows, numColumns, ia, ja, reinterpret_cast<const hmemo::HArray<H>& >( values ), ctx );
        }
        else
        {
            StorageWrapper<Derived, T>::setCSRDataImpl( obj, numRows, numColumns, ia, ja, values, ctx );
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
            StorageWrapper<Derived, T>::buildCSRDataImpl( obj, csrIA, csrJA, csrValues, ctx );
        }
    }

    static void assignImpl(
        Derived* obj,
        const _MatrixStorage& otherStorage )
    {
        if ( otherStorage.getValueType() == common::getScalarType<H>() )
        {
            obj->assignImpl( reinterpret_cast<const MatrixStorage<H>& >( otherStorage ) );
        }
        else
        {
            StorageWrapper<Derived, T>::assignImpl( obj, otherStorage );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
