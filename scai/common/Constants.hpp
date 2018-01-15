/**
 * @file Constants.hpp
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
 * @brief Definition of the Constants One, Zero, Eps for comparisons
 * @author Eric Schricker
 * @date 07.10.2015
 */

#pragma once

#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/exception/InvalidArgumentException.hpp>

#include <cmath>
#include <limits>
#include <iomanip>

namespace scai
{

namespace common
{

/** Own struct for the enum type Constants and its values. */

enum class Constants
{

    /** Enumeration type for constants for which type-specific values are provided */

    ZERO,  //!< stands for value 0
    ONE    //!< stands for value 1

};

/** This method returns the type specific value for each constant */

template<typename ValueType>
inline ValueType getConstant( const Constants& c )
{
    ValueType val( 0 );

    switch ( c )
    {
        case Constants::ONE:
            val = ValueType( 1 );
            break;

        case Constants::ZERO:
            val = ValueType( 0 );
            break;
    }

    return val;
}

/** Comparison against constant ZERO or ONE checks for exact equality */

template<typename ValueType>
inline bool operator==( const ValueType& x, const Constants& c )
{
    if ( Constants::ZERO == c )
    {
        return x == ValueType( 0 );
    }
    else if ( Constants::ONE == c )
    {
        return x == ValueType( 1 );
    }
    else
    {
        SCAI_THROWEXCEPTION( InvalidArgumentException, "Unsupported constant" )
    }
}

// Use template specialization for IndexType as real/imag might not be available

template<>
inline bool operator==( const IndexType& x, const Constants& c )
{
    if ( Constants::ZERO == c )
    {
        return x == 0;
    }
    else
    {
        return x == 1;
    }
}

template<typename ValueType>
bool operator==( const Constants& c, const ValueType& x )
{
    return operator==( x, c );
}

/** Operator not equal also provided for convenience */

template<typename ValueType>
bool operator!=( const ValueType& x, const Constants& c )
{
    return !operator==( x, c );
}

template<typename ValueType>
bool operator!=( const Constants& c, const ValueType& x )
{
    return !operator==( x, c );
}

} /* end namespace common */

} /* end namespace scai */
