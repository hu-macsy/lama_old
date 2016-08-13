/**
 * @file utilskernel/mepr/UtilsWrapper.hpp
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

template<typename TList1, typename TList2> struct UtilsWrapperTT1;
template<typename ValueType, typename TList2> struct UtilsWrapperTT2;

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
    static void setArray( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
    static void gather( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ) {}
    static void scatter( hmemo::HArray<ValueType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
    static void setScalar( hmemo::_HArray&, const ValueType, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
    static void setValImpl( hmemo::_HArray&, const IndexType, const ValueType ) {}
    static ValueType getValImpl( const hmemo::_HArray&, const IndexType )
    {
        return ValueType();
    }
};

template<typename TList> struct UtilsWrapperTT1<common::mepr::NullType, TList>
{
    static void setArray( hmemo::_HArray&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
    static void gather( hmemo::_HArray&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ) {}
    static void scatter( hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
};

template<typename ValueType> struct UtilsWrapperTT2<ValueType, common::mepr::NullType>
{
    static void setArray( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
    static void gather( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ) {}
    static void scatter( hmemo::HArray<ValueType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const reduction::ReductionOp, const hmemo::ContextPtr ) {}
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
    static void setArray( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::setArray( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), op, loc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setArray( target, source, op, loc );
        }
    }

    static void gather( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source,
                        const hmemo::HArray<IndexType>& index, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::gather( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), index, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::gather( target, source, index, prefLoc );
        }
    }

    static void scatter( hmemo::HArray<ValueType>& target,
                         const hmemo::HArray<IndexType>& index,
                         const hmemo::_HArray& source,
                         const reduction::ReductionOp op,
                         const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::scatter( target, index, reinterpret_cast<const hmemo::HArray<H>&>( source ), op, prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::scatter( target, index, source, op, prefLoc );
        }
    }

    static void setScalar( hmemo::_HArray& target, const ValueType value, const reduction::ReductionOp op, const hmemo::ContextPtr ctx )
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

    static void setValImpl( hmemo::_HArray& target, const IndexType index, const ValueType val )
    {
        if ( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setValImpl( reinterpret_cast<hmemo::HArray<H>&>( target ), index, static_cast<H>( val ) );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setValImpl( target, index, val );
        }
    }

    static ValueType getValImpl( const hmemo::_HArray& array, const IndexType index )
    {
        if ( common::getScalarType<H>() == array.getValueType() )
        {
            return static_cast<ValueType>( HArrayUtils::getValImpl( reinterpret_cast<const hmemo::HArray<H>&>( array ), index ) );
        }
        else
        {
            return UtilsWrapperT< ValueType, T>::getValImpl( array, index );
        }
    }
};

template<typename H, typename T, typename TList2>
struct UtilsWrapperTT1<common::mepr::TypeList<H, T>, TList2 >
{
    static void setArray( hmemo::_HArray& target, const hmemo::_HArray& source, const reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperTT2<H, TList2>::setArray( reinterpret_cast<hmemo::HArray<H>&>( target ), source, op, loc );
        }
        else
        {
            UtilsWrapperTT1<T, TList2>::setArray( target, source, op, loc );
        }
    }

    static void gather( hmemo::_HArray& target, const hmemo::_HArray& source, const hmemo::HArray<IndexType>& index, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperTT2<H, TList2>::gather( reinterpret_cast<hmemo::HArray<H>&>( target ), source, index, prefLoc );
        }
        else
        {
            UtilsWrapperTT1<T, TList2>::gather( target, source, index, prefLoc );
        }
    }

    static void scatter( hmemo::_HArray& target, const hmemo::HArray<IndexType>& index, const hmemo::_HArray& source, const reduction::ReductionOp op, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperTT2<H, TList2>::scatter( reinterpret_cast<hmemo::HArray<H>&>( target ), index, source, op, prefLoc );
        }
        else
        {
            UtilsWrapperTT1<T, TList2>::scatter( target, index, source, op, prefLoc );
        }
    }
};

template<typename ValueType, typename H, typename T>
struct UtilsWrapperTT2<ValueType, common::mepr::TypeList<H, T> >
{
    static void setArray( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            HArrayUtils::setArray<ValueType, H>( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), op, loc );
        }
        else
        {
            UtilsWrapperTT2<ValueType, T>::setArray( target, source, op, loc );
        }
    }

    static void gather( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const hmemo::HArray<IndexType>& index, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            HArrayUtils::gather<ValueType, H>( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), index, prefLoc );
        }
        else
        {
            UtilsWrapperTT2<ValueType, T>::gather( target, source, index, prefLoc );
        }
    }

    static void scatter( hmemo::HArray<ValueType>& target, const hmemo::HArray<IndexType>& index, const hmemo::_HArray& source, const reduction::ReductionOp op, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            HArrayUtils::scatter<ValueType, H>( target, index, reinterpret_cast<const hmemo::HArray<H>&>( source ), op, prefLoc );
        }
        else
        {
            UtilsWrapperTT2<ValueType, T>::scatter( target, index, source, op, prefLoc );
        }
    }
};

} /* end namespace mepr */

} /* end namespace utilskernel */

} /* end namespace scai */
