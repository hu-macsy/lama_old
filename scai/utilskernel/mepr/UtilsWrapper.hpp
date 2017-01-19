/**
 * @file utilskernel/mepr/UtilsWrapper.hpp
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
 * @brief Wrapper for templated calls in Utils
 * @author Eric Schricker
 * @date 07.03.2016
 */

#pragma once

// local library
#include <scai/utilskernel/HArrayUtils.hpp>

// scai internal libraries
#include <scai/hmemo/HArray.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/mepr/TypeList.hpp>


namespace scai
{

namespace utilskernel
{

namespace mepr
{

/*
 * Forward declartion
 */
template<typename TList> struct UtilsWrapper;
template<typename ValueType, typename TList> struct UtilsWrapperT;

template<typename TList1, typename TList2> struct UtilsWrapperTT;

/*
 * Termination
 */
template<> struct UtilsWrapper<common::mepr::NullType>
{
    static void setRandom( hmemo::_HArray& array, IndexType, float, hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setRandom not instantiated for " << array.getValueType() )
    }
};

template<typename ValueType> struct UtilsWrapperT<ValueType, common::mepr::NullType>
{
    static void setArray(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArray: source type unsupported, target = " << target << ", source = " << source  )
    }

    static void setArraySection(
        hmemo::_HArray& target,
        const IndexType,
        const IndexType,
        const hmemo::_HArray& source,
        const IndexType,
        const IndexType,
        const IndexType,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArraySection: target type unsupported, target = " << target << ", source = " << source )
    }

    static void gather(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>&,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "gather: source type unsupported, target = " << target << ", source = " << source )
    }

    static void scatter(
        hmemo::HArray<ValueType>& target,
        const hmemo::HArray<IndexType>&,
        const bool,
        const hmemo::_HArray& source,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "scatter: source type unsupported, target = " << target << ", source = " << source )
    }

    static void setScalar( hmemo::_HArray& target, const ValueType, const binary::BinaryOp, const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setScalar: target type unsupported, target = " << target )
    }

    static void setValImpl( hmemo::_HArray& target, const IndexType, const ValueType )
    {
        COMMON_THROWEXCEPTION( "setVal: target type unsupported, target = " << target )
    }

    static ValueType getValImpl( const hmemo::_HArray&, const IndexType )
    {
        return ValueType( 0 );
    }

    static void buildSparse( 
        hmemo::HArray<ValueType>&, 
        hmemo::HArray<IndexType>&,
        const hmemo::_HArray& denseArray, 
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "buildSparse: source type unsupported, array = " << denseArray )
    }
};

template<typename TList> 
struct UtilsWrapperTT<common::mepr::NullType, TList>
{
    static void setArray(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArray: target type unsupported, target = " << target << ", source = " << source )
    }

    static void setArraySection(
        hmemo::_HArray& target,
        const IndexType,
        const IndexType,
        const hmemo::_HArray& source,
        const IndexType,
        const IndexType,
        const IndexType,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArraySection: target type unsupported, target = " << target << ", source = " << source )
    }

    static void gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>&,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "gather: target type unsupported, target = " << target << ", source = " << source )
    }

    static void scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& source,
        const bool,
        const hmemo::_HArray&,
        const binary::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "scatter: target type unsupported, target = " << target << ", source = " << source )
    }

    static void buildSparse( 
        hmemo::_HArray&, 
        hmemo::HArray<IndexType>&,
        const hmemo::_HArray&, 
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "buildSparse: target type unsupported" )
    }
};

/*
 * Step n
 */

template<typename H, typename T>
struct UtilsWrapper<common::mepr::TypeList<H, T> >
{
    static void setRandom( hmemo::_HArray& array, IndexType n, float fillRate, hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() ==  array.getValueType() )
        {
            HArrayUtils::setRandomImpl( reinterpret_cast<hmemo::HArray<H>&>( array ), n, fillRate, loc );
        }
        else
        {
            UtilsWrapper<T>::setRandom( array, n, fillRate, loc );
        }
    }
};

template<typename ValueType, typename H, typename T>
struct UtilsWrapperT< ValueType, common::mepr::TypeList<H, T> >
{
    static void setArray(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const binary::BinaryOp op,
        const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::setArray( target, typedSource, op, loc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setArray( target, source, op, loc );
        }
    }

    static void setArraySection(
        hmemo::HArray<ValueType>& typedTarget, const IndexType targetOffset, const IndexType targetInc,
        const hmemo::_HArray& source, const IndexType sourceOffset, const IndexType sourceInc,
        const IndexType n, const binary::BinaryOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );

            HArrayUtils::setArraySectionImpl<ValueType, H>( typedTarget, targetOffset, targetInc,
                    typedSource, sourceOffset, sourceInc,
                    n, op, loc );
        }
        else
        {
            UtilsWrapperT<ValueType, T>::setArraySection( 
                typedTarget, targetOffset, targetInc,
                source, sourceOffset, sourceInc,
                n, op, loc );
        }
    }

    static void gather(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::gatherImpl( target, typedSource, indexes, op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::gather( target, source, indexes, op, prefLoc );
        }
    }

    static void scatter(
        hmemo::HArray<ValueType>& target,
        const hmemo::HArray<IndexType>& indexes,
        const bool unique,
        const hmemo::_HArray& source,
        const binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::scatterImpl( target, indexes, unique, typedSource, op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::scatter( target, indexes, unique, source, op, prefLoc );
        }
    }

    static void setScalar( hmemo::_HArray& target, const ValueType value, const binary::BinaryOp op, const hmemo::ContextPtr ctx )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setScalar( reinterpret_cast<hmemo::HArray<H>&>( target ), static_cast<H>( value ), op, ctx );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setScalar( target, value, op, ctx );
        }
    }

    static void setValImpl( hmemo::_HArray& target, const IndexType indexes, const ValueType val )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setValImpl( reinterpret_cast<hmemo::HArray<H>&>( target ),
                                     indexes,
                                     static_cast<H>( val ), binary::COPY );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setValImpl( target, indexes, val );
        }
    }

    static ValueType getValImpl( const hmemo::_HArray& array, const IndexType indexes )
    {
        if ( common::getScalarType<H>() == array.getValueType() )
        {
            return static_cast<ValueType>( HArrayUtils::getValImpl( reinterpret_cast<const hmemo::HArray<H>&>( array ), indexes ) );
        }
        else
        {
            return UtilsWrapperT< ValueType, T>::getValImpl( array, indexes );
        }
    }

    static void buildSparse( 
        hmemo::HArray<ValueType>& sparseArray, 
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::_HArray& denseArray, 
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == denseArray.getValueType() )
        {
            const hmemo::HArray<H>& typedDenseArray = reinterpret_cast<const hmemo::HArray<H>&>( denseArray );
            HArrayUtils::buildSparseArrayImpl( sparseArray, sparseIndexes, typedDenseArray, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::buildSparse( sparseArray, sparseIndexes, denseArray, prefLoc );
        }
    }
};

template<typename H, typename T, typename TList2>
struct UtilsWrapperTT<common::mepr::TypeList<H, T>, TList2 >
{
    static void setArray( hmemo::_HArray& target, const hmemo::_HArray& source, const binary::BinaryOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperT<H, TList2>::setArray( reinterpret_cast<hmemo::HArray<H>&>( target ), source, op, loc );
        }
        else
        {
            UtilsWrapperTT<T, TList2>::setArray( target, source, op, loc );
        }
    }

    static void setArraySection( hmemo::_HArray& target, const IndexType targetOffset, const IndexType targetInc,
                                 const hmemo::_HArray& source, const IndexType sourceOffset, const IndexType sourceInc,
                                 const IndexType n, const binary::BinaryOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );

            UtilsWrapperT<H, TList2>::setArraySection( typedTarget, targetOffset, targetInc,
                    source, sourceOffset, sourceInc,
                    n, op, loc );
        }
        else
        {
            UtilsWrapperTT<T, TList2>::setArraySection( target, targetOffset, targetInc,
                    source, sourceOffset, sourceInc,
                    n, op, loc );
        }
    }

    static void gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );
            UtilsWrapperT<H, TList2>::gather( typedTarget, source, indexes, op, prefLoc );
        }
        else
        {
            UtilsWrapperTT<T, TList2>::gather( target, source, indexes, op, prefLoc );
        }
    }

    static void scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& indexes,
        const bool unique,
        const hmemo::_HArray& source,
        const binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );
            UtilsWrapperT<H, TList2>::scatter( typedTarget, indexes, unique, source, op, prefLoc );
        }
        else
        {
            UtilsWrapperTT<T, TList2>::scatter( target, indexes, unique, source, op, prefLoc );
        }
    }

    static void buildSparse(
        hmemo::_HArray& sparseArray,
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::_HArray& denseArray,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == sparseArray.getValueType() )
        {
            hmemo::HArray<H>& typedSparseArray = reinterpret_cast<hmemo::HArray<H>&>( sparseArray );
            UtilsWrapperT<H, TList2>::buildSparse( typedSparseArray, sparseIndexes, denseArray, prefLoc );
        }
        else
        {
            UtilsWrapperTT<T, TList2>::buildSparse( sparseArray, sparseIndexes, denseArray, prefLoc );
        }
    }
};

} /* end namespace mepr */

} /* end namespace utilskernel */

} /* end namespace scai */
