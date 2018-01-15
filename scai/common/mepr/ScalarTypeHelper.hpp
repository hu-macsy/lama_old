/**
 * @file common/mepr/ScalarTypeHelper.hpp
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
 * @brief Conversion from ScalarType to ValueType and other way round
 * @author Eric Schricker
 * @date 10.03.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <typeinfo>

namespace scai
{

namespace common
{

namespace mepr
{

/*
 * Forward declartion
 */

template<typename TList>
struct ScalarTypeHelper;

/*
 * Termination
 */
template<>
struct ScalarTypeHelper<NullType>
{
    static size_t sizeOf( const ScalarType stype )
    {
        COMMON_THROWEXCEPTION( "sizeof: unsupported for " << stype )
        return 0;
    }

    static int precision( const ScalarType )
    {
        // no type T available with TypeTraits<T>::stype
        return 0;
    }

    static bool isComplex( const ScalarType )
    {
        return false;
    }

    static bool isNumeric( const ScalarType )
    {
        return false;
    }

    static bool contains( const ScalarType )
    {
        return false;
    }
};

/*
 * Step n
 */
template<typename H, typename T>
struct ScalarTypeHelper< TypeList<H, T> >
{
    static size_t sizeOf( const ScalarType s )
    {
        if ( s == common::getScalarType<H>() )
        {
            return sizeof( H );
        }
        else
        {
            return ScalarTypeHelper<T>::sizeOf( s );
        }
    }

    static int precision( const ScalarType stype )
    {
        if ( stype == TypeTraits<H>::stype )
        {
            return TypeTraits<H>::precision();
        }
        else
        {
            return ScalarTypeHelper<T>::precision( stype );
        }
    }

    static bool isComplex( const ScalarType stype )
    {
        if ( stype == TypeTraits<H>::stype )
        {
            // a type is complex if its RealType is not the same

            typedef typename TypeTraits<H>::RealType RealType;

            return typeid( H ) != typeid( RealType );
        }
        else
        {
            return ScalarTypeHelper<T>::isComplex( stype );
        }
    }

    static bool isNumeric( const ScalarType stype )
    {
        if ( stype == TypeTraits<H>::stype )
        {
            // a type is numeric if its small value is not ZERO

            return TypeTraits<H>::small() != H( 0 );
        }
        else
        {
            return ScalarTypeHelper<T>::isNumeric( stype );
        }
    }

    static bool contains( const ScalarType stype )
    {
        if ( stype == TypeTraits<H>::stype )
        {
            return true;
        }
        else
        {
            return ScalarTypeHelper<T>::contains( stype );
        }
    }

};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
