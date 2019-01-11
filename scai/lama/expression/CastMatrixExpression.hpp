/**
 * @file CastMatrixExpression.hpp
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
 * @brief Functions to build symbolic expressions cast<TargetValueType>( Matrix<SourceValueType> )
 * @author Thomas Brandes
 * @date 22.12.2017
 */
#pragma once

namespace scai
{

namespace lama
{

template<typename ValueType>
class Matrix;

/* ============================================================================ */
/*    CastMatrixExpression                                                      */
/* ============================================================================ */

/** Class that is used for symbolic expressions to cast matrices from one
 *  type to another type.
 *
 *  @tparam TargetType is the value type after conversion
 *  @tparam SourceType is the value type before conversion
 */
template<typename TargetType, typename SourceType>
class CastMatrixExpression
{
public:

    CastMatrixExpression( const Matrix<SourceType>& arg ) :  mArg( arg )
    {
    }

    inline const Matrix<SourceType>& getArg() const
    {
        return mArg;
    }

private:

    const Matrix<SourceType>& mArg;
};

    
/* ============================================================================ */
/*    functions tto create CastMatrixExpression's                               */
/* ============================================================================ */

/**
 *  @brief Build symbolic expression to cast matrix from one type (SourceType) to another (TargetType).
 *
 * \code
 *    CSRSparseMatrix<float> fV;
 *    DenseMatrix<double> dV;
 *    fV = dV;                           // illegal 
 *    fV = cast<float, double>( dv );    // convert double -> float
 *    dV = cast<double, float>( fv );    // convert float -> double
 *
 *    dV = cast<double>( fv );    // okay, 2nd arg float can be deduced
 *    fV = cast<float>( fv );     // okay, 2nd arg double can be deduced
 * \endcode
 */
template<typename TargetType, typename SourceType>
CastMatrixExpression<TargetType, SourceType> cast( const Matrix<SourceType>& v )
{
    return CastMatrixExpression<TargetType, SourceType>( v );
}

} /* end namespace lama */

} /* end namespace scai */

