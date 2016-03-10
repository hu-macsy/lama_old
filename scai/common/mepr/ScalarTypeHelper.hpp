/*
 * ScalarTypeHelper.hpp
 *
 *  Created on: Mar 10, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai {

namespace common {

namespace mepr {

template<typename TList>
struct ScalarTypeHelper;

template<>
struct ScalarTypeHelper<NullType>
{
    static long sizeOf( const scalar::ScalarType& )
    {
        return 0;
    }

    static scalar::ScalarType getBySize( const long )
    {
        return scalar::UNKNOWN;
    }
};

template<typename H, typename T>
struct ScalarTypeHelper< TypeList<H,T> >
{
    static long sizeOf( const scalar::ScalarType& s )
    {
        if( s == common::getScalarType<H>() )
        {
            return sizeof( H );
        }
        else
        {
            return ScalarTypeHelper<T>::sizeOf( s );
        }
    }

    static scalar::ScalarType getBySize( const long size )
    {
        if( sizeof( H ) == size )
        {
            return TypeTraits<H>::stype;
        }
        else
        {
            return ScalarTypeHelper< T >::getBySize( size );
        }
    }
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
