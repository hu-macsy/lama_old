/**
 * @file CompareOp.hpp
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
 * @brief Enumeration type for the different comparison operators
 * @author Thomas Brandes
 * @date 15.01.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/throw.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
{

/** 
 * Enumeration type for binary compare operators used in array operations.
 */
enum class CompareOp 
{
    LT,            //!< for less than
    LE,            //!< for less equal
    GE,            //!< for greater equal
    GT,            //!< for greater than
    EQ,            //!< for equality
    NE,            //!< for non equal
    MAX_COMPARE_OP //!< for internal use only
};

/** 
 * @brief function that compares two values corresponding to the specified operator.
 *
 * Note: this method should not be used for complex values (at runtime) but it is
 *       provided in such a way that it can be compiled in other code to avoid template specializations.
 */
template <typename ValueType>
CUDA_CALLABLE_MEMBER
inline bool compare( const ValueType& x1, const CompareOp op, const ValueType& x2 )
{
    // just imagine here an assert that ValueType is not complex
    // but you don't see it here as exceptions might not work on device code

    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    RealType v1 = x1;
    RealType v2 = x2;

    switch ( op )
    {
        case CompareOp::LT:
            return v1 < v2;
        case CompareOp::LE:
            return v1 <= v2;
        case CompareOp::GE:
            return v1 >= v2;
        case CompareOp::GT:
            return v1 > v2;
        case CompareOp::EQ:
            return v1 == v2;
        case CompareOp::NE:
            return v1 != v2;
        default:
            return false;
    }
}

/*
 * Output of CompareOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const CompareOp& op )
{
    switch ( op )
    {
        case CompareOp::LE:
            stream << "LE";
            break;

        case CompareOp::LT:
            stream << "LT";
            break;

        case CompareOp::GE:
            stream << "GE";
            break;

        case CompareOp::GT:
            stream << "GT";
            break;

        case CompareOp::EQ:
            stream << "EQ";
            break;

        case CompareOp::NE:
            stream << "NE";
            break;

        default:
            stream << "<unknown_compare_op>";
            break;
    }

    return stream;
}

} /* end namespace common */

} /* end namespace scai */
