/**
 * @file CastVectorExpression.hpp
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
 * @brief Functions to build symbolic expressions for conversion between vectors.
 * @author Thomas Brandes
 * @date 22.12.2017
 */
#pragma once

namespace scai
{

namespace lama
{

template<typename ValueType>
class Vector;

/* ============================================================================ */
/*    CastVectorExpression                                                      */
/* ============================================================================ */

/** Class that is used for symbolic expressions to cast vectors from one
 *  type to other type.
 *
 *  @tparam TargetType is the value type after conversion
 *  @tparam SourceType is the value type before conversion
 */
template<typename TargetType, typename SourceType>
class CastVectorExpression
{
public:

    CastVectorExpression( const Vector<SourceType>& arg ) :  mArg( arg )
    {
    }

    inline const Vector<SourceType>& getArg() const
    {
        return mArg;
    }

private:

    const Vector<SourceType>& mArg;
};

    
/* ============================================================================ */
/*    functions tto create CastVectorExpression's                               */
/* ============================================================================ */

/**
 *  Free function to cast vector from one type (SourceType) to another (TargetType).
 *
 * \code
 *    DenseVector<float> fV;
 *    DenseVector<double> dV;
 *    fV = dV;                           // illegal 
 *    fV = cast<float, double>( dv );    // convert double -> float
 *    dV = cast<double, float>( fv );    // convert float -> double
 *
 *    dV = cast<double>( fv );    // okay, 2nd arg float can be deduced
 *    fV = cast<float>( fv );     // okay, 2nd arg double can be deduced
 * \endcode
 */
template<typename TargetType, typename SourceType>
CastVectorExpression<TargetType, SourceType> cast( const Vector<SourceType>& v )
{
    return CastVectorExpression<TargetType, SourceType>( v );
}

} /* end namespace lama */

} /* end namespace scai */

