/**
 * @file ComplexVectorExpression.hpp
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
 * @brief Functions to build symbolic expressions for compound and select complex vectors.
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

/* ============================================================================= */
/*    ComplexPartVectorExpression, e.g. imag( DenseVector<ComplexFloat> x ) */
/* ============================================================================= */

/** 
 *  Expression to select real or imaginary part of a vector
 *
 *  @tparam ValueType is the type of the (complex) vector that is selected
 *  @tparam kind      is either REAL or IMAG correspodinng to the selected part
 */
template<typename ValueType, common::ComplexPart kind>
class ComplexPartVectorExpression
{
public:

    ComplexPartVectorExpression( const Vector<ValueType>& arg ) :  mArg( arg )
    {
    }

    inline const Vector<ValueType>& getArg() const
    {
        return mArg;
    }

private:

    const Vector<ValueType>& mArg;   // reference to the function argument
};

template<typename ValueType>
ComplexPartVectorExpression<ValueType, common::ComplexPart::REAL> real( const Vector<ValueType>& v )
{
    return ComplexPartVectorExpression<ValueType, common::ComplexPart::REAL>( v );
}

template<typename ValueType>
ComplexPartVectorExpression<ValueType, common::ComplexPart::IMAG> imag( const Vector<ValueType>& v )
{
    return ComplexPartVectorExpression<ValueType, common::ComplexPart::IMAG>( v );
}

/** A complex vector expression is very close to a binary expression that stands for x + y * i.
 *  In contrary to other binary operations we have a new result type so it is handled here on its own
 */
template<typename ValueType>
class ComplexBuildVectorExpression
{
public:

    ComplexBuildVectorExpression( const Vector<ValueType>& x, const Vector<ValueType>& y ) :  

        mReal( x ),
        mImag( y )
    {
    }

    inline const Vector<ValueType>& getRealArg() const
    {
        return mReal;
    }

    inline const Vector<ValueType>& getImagArg() const
    {
        return mImag;
    }

private:

    const Vector<ValueType>& mReal;
    const Vector<ValueType>& mImag;
};

/** Method to build symbolic expression for complex( x, y ), stands for elementwise x + i * y  */

template<typename ValueType>
ComplexBuildVectorExpression<ValueType> complex( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return ComplexBuildVectorExpression<ValueType>( x, y );
}

} /* end namespace lama */

} /* end namespace scai */

