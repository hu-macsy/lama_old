/**
 * @file utilskernel/mepr/UtilsWrapper.hpp
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


namespace scai {

namespace utilskernel {

namespace mepr {

/*
 * Forward declartion
 */
template<typename ValueType, typename TList> struct UtilsWrapperT;

template<typename TList1, typename TList2> struct UtilsWrapperTT1;
template<typename ValueType, typename TList2> struct UtilsWrapperTT2;

/*
 * Termination
 */
template<typename ValueType> struct UtilsWrapperT<ValueType,common::mepr::NullType>
{
    static void setArray( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void gather( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ){}
    static void scatter( hmemo::HArray<ValueType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const hmemo::ContextPtr ){}
    static void setScalar( hmemo::_HArray&, const ValueType, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void setValImpl( hmemo::_HArray&, const IndexType, const ValueType ){}
    static ValueType getValImpl( const hmemo::_HArray&, const IndexType ){ return ValueType(); }
};

template<typename TList> struct UtilsWrapperTT1<common::mepr::NullType, TList>
{
    static void setArray( hmemo::_HArray&, const hmemo::_HArray&, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void gather( hmemo::_HArray&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ){}
    static void scatter( hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const hmemo::ContextPtr ){}
};

template<typename ValueType> struct UtilsWrapperTT2<ValueType, common::mepr::NullType>
{
    static void setArray( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void gather( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const hmemo::HArray<IndexType>&, const hmemo::ContextPtr ){}
    static void scatter( hmemo::HArray<ValueType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&, const hmemo::ContextPtr ){}
};

/*
 * Step n
 */
template<typename ValueType, typename H, typename T>
struct UtilsWrapperT< ValueType, common::mepr::TypeList<H, T> >
{
    static void setArray( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if( common::getScalarType<H>() ==  source.getValueType() )
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
        if( common::getScalarType<H>() ==  source.getValueType() )
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
                         const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::scatter( target, index, reinterpret_cast<const hmemo::HArray<H>&>( source ), prefLoc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::scatter( target, index, source, prefLoc );
        }
    }

    static void setScalar( hmemo::_HArray& target, const ValueType value, const common::reduction::ReductionOp op, const hmemo::ContextPtr ctx )
    {
        if( common::getScalarType<H>() ==  target.getValueType() )
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
        if( common::getScalarType<H>() ==  target.getValueType() )
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
        if( common::getScalarType<H>() == array.getValueType() )
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
    static void setArray( hmemo::_HArray& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if( common::getScalarType<H>() == target.getValueType() )
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
        if( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperTT2<H, TList2>::gather( reinterpret_cast<hmemo::HArray<H>&>( target ), source, index, prefLoc );
        }
        else
        {
            UtilsWrapperTT1<T, TList2>::gather( target, source, index, prefLoc );
        }
    }

    static void scatter( hmemo::_HArray& target, const hmemo::HArray<IndexType>& index, const hmemo::_HArray& source, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperTT2<H, TList2>::scatter( reinterpret_cast<hmemo::HArray<H>&>( target ), index, source, prefLoc );
        }
        else
        {
            UtilsWrapperTT1<T, TList2>::scatter( target, index, source, prefLoc );
        }
    }
};

template<typename ValueType, typename H, typename T>
struct UtilsWrapperTT2<ValueType, common::mepr::TypeList<H, T> >
{
    static void setArray( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if( common::getScalarType<H>() == source.getValueType() )
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
        if( common::getScalarType<H>() == source.getValueType() )
        {
            HArrayUtils::gather<ValueType, H>( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), index, prefLoc );
        }
        else
        {
            UtilsWrapperTT2<ValueType, T>::gather( target, source, index, prefLoc );
        }
    }

    static void scatter( hmemo::HArray<ValueType>& target, const hmemo::HArray<IndexType>& index, const hmemo::_HArray& source, const hmemo::ContextPtr prefLoc )
    {
        if ( common::getScalarType<H>() == source.getValueType() )
        {
            HArrayUtils::scatter<ValueType, H>( target, index, reinterpret_cast<const hmemo::HArray<H>&>( source ), prefLoc );
        }
        else
        {
            UtilsWrapperTT2<ValueType, T>::scatter( target, index, source, prefLoc );
        }
    }
};

} /* end namespace mepr */

} /* end namespace utilskernel */

} /* end namespace scai */
