/**
 * @file utilskernel/mepr/UtilsWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

/* ------------------------------------------------------------------------- */
/*  Forward declaration of used structures                                   */
/* ------------------------------------------------------------------------- */

/** This struct is used as metaprogramming to call a routine for each type of a list
 *
 *  @tparam TList is a type list containing all types that will be considered.
 */
template<typename TList> struct UtilsWrapper;

/** This struct is used as metaprogramming to call a template routine for each type of a list
 *
 *  @tparam ValueType is an addition type used for the methods that will be called
 *  @tparam TList is a type list containing all types that will be considered.
 */
template<typename ValueType, typename TList> struct UtilsWrapperT;

/** This struct is as metaprogramm to call a template routine for each pair of types of two type lists.
 *
 *  @tparam TList1 is the first type list
 *  @tparam TList2 is the seond type list
 *
 *  This routine traverses at first TList1 and calls for each type T UtilsWrapperT<T, TList2>.
 */
template<typename TList1, typename TList2> struct UtilsWrapperTT;

/* ------------------------------------------------------------------------- */
/*  Termination                                                              */
/* ------------------------------------------------------------------------- */

template<> struct UtilsWrapper<common::mepr::NullType>
{
    static void fillRandom( hmemo::_HArray& array, IndexType, float, hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "fillRandom not instantiated for " << array.getValueType() )
    }
};

template<typename ValueType> 
struct UtilsWrapperT<ValueType, common::mepr::NullType>
{
    static void setArray(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArray: source type " << source.getValueType() << " unsupported, target = " 
                               << target << ", source = " << source  )
    }

    static void setArraySection(
        hmemo::_HArray& target,
        const IndexType,
        const IndexType,
        const hmemo::_HArray& source,
        const IndexType,
        const IndexType,
        const IndexType,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArraySection: target type unsupported, target = " << target << ", source = " << source )
    }

    static void gather(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>&,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "gather: source type unsupported, target = " << target << ", source = " << source )
    }

    static void sparseGather(
        hmemo::_HArray& target,
        const ValueType,
        const hmemo::HArray<ValueType>&,
        const hmemo::HArray<IndexType>&,
        const hmemo::HArray<IndexType>&,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "sparseGather: target type " << target.getValueType() << " unsupported" )
    }

    static void scatter(
        hmemo::HArray<ValueType>& target,
        const hmemo::HArray<IndexType>&,
        const bool,
        const hmemo::_HArray& source,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "scatter: no instantiation for type " << source.getValueType() << " of source array"
                               << ", target = " << target )
    }

    static void setScalar( hmemo::_HArray& target, const ValueType, const common::BinaryOp, const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setScalar: target type unsupported, target = " << target )
    }

    static void setSameValue( hmemo::_HArray& target, const IndexType, const ValueType, const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setSameValue: target type unsupported, target = " << target )
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
        const common::BinaryOp,
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
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "setArraySection: target type unsupported, target = " << target << ", source = " << source )
    }

    static void gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>&,
        const common::BinaryOp,
        const hmemo::ContextPtr )
    {
        COMMON_THROWEXCEPTION( "gather: target type unsupported, target = " << target << ", source = " << source )
    }

    static void scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& source,
        const bool,
        const hmemo::_HArray&,
        const common::BinaryOp,
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

/* ------------------------------------------------------------------------- */
/*  Termination                                                              */
/* ------------------------------------------------------------------------- */

template<typename THead, typename TList>
struct UtilsWrapper<common::mepr::TypeList<THead, TList> >
{
    static void fillRandom( hmemo::_HArray& array, IndexType bound, float fillRate, hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<THead>() ==  array.getValueType() )
        {
            // now type of array is known, so we can make a safe cast
            hmemo::HArray<THead>& typedArray = reinterpret_cast<hmemo::HArray<THead>&>( array );
            HArrayUtils::fillRandom( typedArray, bound, fillRate, loc );
        }
        else
        {
            UtilsWrapper<TList>::fillRandom( array, bound, fillRate, loc );
        }
    }
};

template<typename ValueType, typename H, typename Tail>
struct UtilsWrapperT< ValueType, common::mepr::TypeList<H, Tail> >
{
    static void setArray(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const common::BinaryOp op,
        const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::setArray( target, typedSource, op, loc );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::setArray( target, source, op, loc );
        }
    }

    static void setArraySection(
        hmemo::HArray<ValueType>& typedTarget, const IndexType targetOffset, const IndexType targetInc,
        const hmemo::_HArray& source, const IndexType sourceOffset, const IndexType sourceInc,
        const IndexType n, const common::BinaryOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );

            HArrayUtils::setArraySection<ValueType, H>( typedTarget, targetOffset, targetInc,
                    typedSource, sourceOffset, sourceInc,
                    n, op, loc );
        }
        else
        {
            UtilsWrapperT<ValueType, Tail>::setArraySection( 
                typedTarget, targetOffset, targetInc,
                source, sourceOffset, sourceInc,
                n, op, loc );
        }
    }

    static void gather(
        hmemo::HArray<ValueType>& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::gather( target, typedSource, indexes, op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::gather( target, source, indexes, op, prefLoc );
        }
    }

    static void sparseGather(
        hmemo::_HArray& target,
        const ValueType sourceZ,
        const hmemo::HArray<ValueType>& sourceVA,
        const hmemo::HArray<IndexType>& sourceIA,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );
            HArrayUtils::sparseGather( typedTarget, sourceZ, sourceVA, sourceIA, indexes, op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::sparseGather( target, sourceZ, sourceVA, sourceIA, indexes, op, prefLoc );
        }
    }

    static void scatter(
        hmemo::HArray<ValueType>& target,
        const hmemo::HArray<IndexType>& indexes,
        const bool unique,
        const hmemo::_HArray& source,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            const hmemo::HArray<H>& typedSource = reinterpret_cast<const hmemo::HArray<H>&>( source );
            HArrayUtils::scatter( target, indexes, unique, typedSource, op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::scatter( target, indexes, unique, source, op, prefLoc );
        }
    }

    static void setScalar( hmemo::_HArray& target, const ValueType value, const common::BinaryOp op, const hmemo::ContextPtr ctx )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setScalar( reinterpret_cast<hmemo::HArray<H>&>( target ), static_cast<H>( value ), op, ctx );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::setScalar( target, value, op, ctx );
        }
    }

    static void setSameValue( hmemo::_HArray& target, const IndexType size, const ValueType val, const hmemo::ContextPtr ctx )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setSameValue( static_cast<hmemo::HArray<H>&>( target ), size, static_cast<H>( val ), ctx );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::setSameValue( target, size, val, ctx );
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
            H zero = 0;
            HArrayUtils::buildSparseArray( sparseArray, sparseIndexes, typedDenseArray, zero, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, Tail >::buildSparse( sparseArray, sparseIndexes, denseArray, prefLoc );
        }
    }
};

template<typename H, typename Tail, typename TList2>
struct UtilsWrapperTT<common::mepr::TypeList<H, Tail>, TList2 >
{
    static void setArray( hmemo::_HArray& target, const hmemo::_HArray& source, const common::BinaryOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperT<H, TList2>::setArray( reinterpret_cast<hmemo::HArray<H>&>( target ), source, op, loc );
        }
        else
        {
            UtilsWrapperTT<Tail, TList2>::setArray( target, source, op, loc );
        }
    }

    static void setArraySection( hmemo::_HArray& target, const IndexType targetOffset, const IndexType targetInc,
                                 const hmemo::_HArray& source, const IndexType sourceOffset, const IndexType sourceInc,
                                 const IndexType n, const common::BinaryOp op, const hmemo::ContextPtr loc )
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
            UtilsWrapperTT<Tail, TList2>::setArraySection( target, targetOffset, targetInc,
                    source, sourceOffset, sourceInc,
                    n, op, loc );
        }
    }

    static void gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );
            UtilsWrapperT<H, TList2>::gather( typedTarget, source, indexes, op, prefLoc );
        }
        else
        {
            UtilsWrapperTT<Tail, TList2>::gather( target, source, indexes, op, prefLoc );
        }
    }

    static void scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& indexes,
        const bool unique,
        const hmemo::_HArray& source,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            hmemo::HArray<H>& typedTarget = reinterpret_cast<hmemo::HArray<H>&>( target );
            UtilsWrapperT<H, TList2>::scatter( typedTarget, indexes, unique, source, op, prefLoc );
        }
        else
        {
            UtilsWrapperTT<Tail, TList2>::scatter( target, indexes, unique, source, op, prefLoc );
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
            UtilsWrapperTT<Tail, TList2>::buildSparse( sparseArray, sparseIndexes, denseArray, prefLoc );
        }
    }
};

} /* end namespace mepr */

} /* end namespace utilskernel */

} /* end namespace scai */
