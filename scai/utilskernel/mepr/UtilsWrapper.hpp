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

template<typename TList> struct UtilsWrapper;

/*
 * Termination
 */
template<typename ValueType> struct UtilsWrapperT<ValueType,common::mepr::NullType>
{
    static void setImpl( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void gatherImpl( hmemo::HArray<ValueType>&, const hmemo::_HArray&, const hmemo::HArray<IndexType>& ){}
    static void setScalarImpl( hmemo::_HArray&, const ValueType, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void setValImpl( hmemo::_HArray&, const IndexType, const ValueType ){}
    static ValueType getValImpl( const hmemo::_HArray&, const IndexType ){ return ValueType(); }
};

template<> struct UtilsWrapper<common::mepr::NullType>
{
    static void setImpl( hmemo::_HArray&, const hmemo::_HArray&, const common::reduction::ReductionOp, const hmemo::ContextPtr ){}
    static void gatherImpl( hmemo::_HArray&, const hmemo::_HArray&, const hmemo::HArray<IndexType>& ){}
};


/*
 * Step n
 */
template<typename ValueType, typename H, typename T>
struct UtilsWrapperT< ValueType, common::mepr::TypeList<H,T> >
{
    static void setImpl( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::setImpl( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), op, loc );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setImpl( target, source, op, loc );
        }
    }

    static void gatherImpl( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const hmemo::HArray<IndexType>& index )
    {
        if( common::getScalarType<H>() ==  source.getValueType() )
        {
            HArrayUtils::gatherImpl( target, reinterpret_cast<const hmemo::HArray<H>&>( source ), index );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::gatherImpl( target, source, index );
        }
    }

    static void setScalarImpl( hmemo::_HArray& target, const ValueType value, const common::reduction::ReductionOp op, const hmemo::ContextPtr ctx )
    {
        if( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setScalarImpl( reinterpret_cast<hmemo::HArray<H>&>( target ), value, op, ctx );
        }
        else
        {
            UtilsWrapperT< ValueType, T >::setScalarImpl( target, value, op, ctx );
        }
    }

    static void setValImpl( hmemo::_HArray& target, const IndexType index, const ValueType val )
    {
        if( common::getScalarType<H>() ==  target.getValueType() )
        {
            HArrayUtils::setValImpl( reinterpret_cast<hmemo::HArray<H>&>( target ), index, val );
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
            return HArrayUtils::getValImpl< ValueType, H >( reinterpret_cast<const hmemo::HArray<H>&>( array ), index );
        }
        else
        {
            return UtilsWrapperT< ValueType, T>::getValImpl( array, index );
        }
    }
};

template<typename H, typename T>
struct UtilsWrapper<common::mepr::TypeList<H,T> >
{
    static void setImpl( hmemo::_HArray& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
    {
        if( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperT<H, ARITHMETIC_ARRAY_HOST_LIST>::setImpl( reinterpret_cast<hmemo::HArray<H>&>( target ), source, op, loc );
        }
        else
        {
            UtilsWrapper<T>::setImpl( target, source, op, loc );
        }
    }

    static void gatherImpl( hmemo::_HArray& target, const hmemo::_HArray& source, const hmemo::HArray<IndexType> index )
    {
        if( common::getScalarType<H>() == target.getValueType() )
        {
            UtilsWrapperT<H, ARITHMETIC_ARRAY_HOST_LIST>::gatherImpl( reinterpret_cast<hmemo::HArray<H>&>( target ), source, index );
        }
        else
        {
            UtilsWrapper<T>::gatherImpl( target, source, index );
        }
    }
};

} /* end namespace mepr */

} /* end namespace utilskernel */

} /* end namespace scai */
