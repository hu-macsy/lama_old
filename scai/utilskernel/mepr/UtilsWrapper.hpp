/*
 * UtilsWrapper.hpp
 *
 *  Created on: Mar 7, 2016
 *      Author: eschricker
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
