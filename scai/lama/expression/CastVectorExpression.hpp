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
 * @brief Functions to build symbolic expressions cast<TargetValueType>( Vector<SourceValueType> )
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
 *  type to ther type.
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

/** 
 *  Expression to select real or imaginary part of a vector
 *  
 *  Note: the input argument must be a complex vector
 */
template<typename RealValueType, common::ComplexSelection kind>
class SelectVectorExpression
{
public:

    SelectVectorExpression( const Vector<common::Complex<RealValueType> >& arg ) :  mSelection( kind ), mArg( arg )
    {
    }

    inline const Vector<common::Complex<RealValueType>>& getArg() const
    {
        return mArg;
    }

private:

    const common::ComplexSelection mSelection;
    const Vector<common::Complex<RealValueType> >& mArg;
};

template<typename RealValueType>
SelectVectorExpression<RealValueType, common::ComplexSelection::REAL> real( const Vector<common::Complex<RealValueType> >& v )
{
    return SelectVectorExpression<RealValueType, common::ComplexSelection::REAL>( v );
}

template<typename RealValueType>
SelectVectorExpression<RealValueType, common::ComplexSelection::IMAG> imag( const Vector<common::Complex<RealValueType> >& v )
{
    return SelectVectorExpression<RealValueType, common::ComplexSelection::IMAG>( v );
}

/** A complex vector expression is very close to a binary expression that stands for x + y * i.
 *  In contrary to other binary operations we have a new result type so it is handled here differently.
 */
template<typename RealValueType>
class ComplexVectorExpression
{
public:

    ComplexVectorExpression( const Vector<RealValueType>& x, const Vector<RealValueType>& y ) :  

        mReal( x ),
        mImag( y )
    {
    }

    inline const Vector<RealValueType>& getRealArg() const
    {
        return mReal;
    }

    inline const Vector<RealValueType>& getImagArg() const
    {
        return mImag;
    }

private:

    const Vector<RealValueType>& mReal;
    const Vector<RealValueType>& mImag;
};

/* cmplx( x, y ) stands for elementwise x + i * y  */

template<typename RealValueType>
ComplexVectorExpression<RealValueType> cmplx( const Vector<RealValueType>& x, const Vector<RealValueType>& y )
{
    return ComplexVectorExpression<RealValueType>( x, y );
}

} /* end namespace lama */

} /* end namespace scai */

