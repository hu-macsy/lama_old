/**
 * @file common/mepr/ScalarTypeHelper.hpp
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

namespace scai {

namespace common {

namespace mepr {

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
    static long sizeOf( const scalar::ScalarType& )
    {
        return 0;
    }

    static scalar::ScalarType getBySize( const long )
    {
        return scalar::UNKNOWN;
    }
};

/*
 * Step n
 */
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
