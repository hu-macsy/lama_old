/**
 * @file ComplexMatrixExpression.hpp
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
 * @brief Functions to build symbolic expressions for compound and select parts of complex matrices.
 * @author Thomas Brandes
 * @date 22.12.2017
 */
#pragma once

#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class Matrix;

/* ============================================================================= */
/*    ComplexPartMatrixExpression, e.g. imag( DenseMatrix<ComplexFloat> x ) */
/* ============================================================================= */

/** 
 *  Expression to select real or imaginary part of a matrix
 *
 *  @tparam ValueType is the type of the (complex) matrix that is selected
 *  @tparam kind      is either REAL or IMAG correspodinng to the selected part
 */
template<typename ValueType, common::ComplexPart kind>
class ComplexPartMatrixExpression
{
public:

    ComplexPartMatrixExpression( const Matrix<ValueType>& arg ) :  mArg( arg )
    {
    }

    inline const Matrix<ValueType>& getArg() const
    {
        return mArg;
    }

private:

    const Matrix<ValueType>& mArg;   // reference to the function argument
};

/** 
 *  @brief Build a symbolic expression for real( complexMatrix ) 
 *
 *  Like all other symbolic expressions in LAMA this expression is resolved by an assignment
 *  or constructor of a matrix.
 */
template<typename ValueType>
ComplexPartMatrixExpression<ValueType, common::ComplexPart::REAL> real( const Matrix<ValueType>& v )
{
    return ComplexPartMatrixExpression<ValueType, common::ComplexPart::REAL>( v );
}

template<typename ValueType>
ComplexPartMatrixExpression<ValueType, common::ComplexPart::IMAG> imag( const Matrix<ValueType>& v )
{
    return ComplexPartMatrixExpression<ValueType, common::ComplexPart::IMAG>( v );
}

/** 
 *  @brief Symbolix expression class for compounding two real matrices to a complex matrix. 
 *
 *  A complex matrix expression is very close to a binary expression that stands for matrix1 + matrix2 * i.
 *  In contrary to other binary operations we have a new result type so it is handled here on its own
 */
template<typename ValueType>
class ComplexBuildMatrixExpression
{
public:

    ComplexBuildMatrixExpression( const Matrix<ValueType>& x, const Matrix<ValueType>& y ) :  

        mReal( x ),
        mImag( y )
    {
    }

    inline const Matrix<ValueType>& getRealArg() const
    {
        return mReal;
    }

    inline const Matrix<ValueType>& getImagArg() const
    {
        return mImag;
    }

private:

    const Matrix<ValueType>& mReal;
    const Matrix<ValueType>& mImag;
};

/** Method to build symbolic expression for complex( x, y ), stands for elementwise x + i * y  */

template<typename ValueType>
ComplexBuildMatrixExpression<ValueType> complex( const Matrix<ValueType>& x, const Matrix<ValueType>& y )
{
    return ComplexBuildMatrixExpression<ValueType>( x, y );
}

} /* end namespace lama */

} /* end namespace scai */

