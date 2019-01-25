/**
 * @file ExpressionTraits.hpp
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
 * @brief Structure to deduce the value type of matrix/vector expressions
 * @author Thomas Brandes
 * @date 25.01.2019
 */
#pragma once

#include <scai/lama/expression/UnaryVectorExpression.hpp>
#include <scai/lama/expression/CastVectorExpression.hpp>
#include <scai/lama/expression/ComplexVectorExpression.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

namespace scai
{

namespace lama
{

/**
 *  @brief Template struct that allows to deduce the ValueType of an expression
 *
 *  This is the general/default definition.
 */
template<typename expression>
struct COMMON_DLL_IMPORTEXPORT ExpressionTraits
{
};

/**
 *  @brief Deduce the ValueType of a vector.
 */
template<typename ValueType>
struct ExpressionTraits<Vector<ValueType>>
{
    typedef ValueType ExpType;
};

/**
 *  @brief Deduce the ValueType of a matrix.
 */
template<typename ValueType>
struct ExpressionTraits<Matrix<ValueType>>
{
    typedef ValueType ExpType;
};

/**
 *  @brief Deduce the value type of a binary expressions, do not take scalar 
 */
template<typename exp, common::BinaryOp op>
struct ExpressionTraits< Expression<exp, intern::Scalar, op> >
{
    typedef typename ExpressionTraits<exp>::ExpType ExpType;
};

/**
 *  @brief Deduce the value type of a binary expressions, do not take scalar 
 */
template<typename exp, common::BinaryOp op>
struct ExpressionTraits< Expression<intern::Scalar, exp, op> >
{
    typedef typename ExpressionTraits<exp>::ExpType ExpType;
};

/**
 *  @brief Deduce the value type of a binary expressions, just take one
 */
template<typename exp1, typename exp2, common::BinaryOp op>
struct ExpressionTraits< Expression<exp1, exp2, op> >
{
    typedef typename ExpressionTraits<exp2>::ExpType ExpType;
};

/**
 *  @brief Type of a matrix operation (transpose, conj), do not change the type.
 */
template<typename ValueType>
struct ExpressionTraits< OpMatrix<ValueType> >
{
    typedef ValueType ExpType;
};

/**
 *  @brief Selection of a complex part gives the real part, also as type.
 */
template<typename ValueType, common::ComplexPart kind>
struct ExpressionTraits< ComplexPartVectorExpression<ValueType, kind> >
{
    typedef RealType<ValueType> ExpType;
};

/**
 *  @brief Building a complex vector gives the complex type
 */
template<typename ValueType>
struct ExpressionTraits< ComplexBuildVectorExpression<ValueType> >
{
    typedef common::Complex<RealType<ValueType>> ExpType;
};

/**
 *  @brief Casting of a vector
 */
template<typename TargetType, typename SourceType>
struct ExpressionTraits< CastVectorExpression<TargetType, SourceType> >
{
    typedef TargetType ExpType;
};

} /* end namespace lama */

} /* end namespace scai */
