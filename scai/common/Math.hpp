/**
 * @file Math.hpp
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
 * @brief Definition of struct that contains math routines overloaded for all SCAI types
 * @author Eric Schricker
 * @date 26.01.2016
 */

#pragma once

#include <scai/common/Math.hpp>

#ifdef SCAI_COMPLEX_SUPPORTED
#include <scai/common/Complex.hpp>
#endif

namespace scai
{

namespace common
{

/** 
 *  Math = MathReal + same routines for complex numbers
 */

struct Math : public MathReal 
{

#ifdef SCAI_COMPLEX_SUPPORTED

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sqrt( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sqrt( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sqrt( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double abs( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> conj( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> conj( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> conj( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> exp( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> exp( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> exp( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> pow( const Complex<float>& base, const Complex<float>& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> pow( const Complex<double>& base, const Complex<double>& exponent );

    static inline MIC_CALLABLE_MEMBER Complex<long double> pow( const Complex<long double>& base, const Complex<long double>& exponent );
 
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> mod( const Complex<float>& x, const Complex<float>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> mod( const Complex<double>& x, const Complex<double>& y );

    static inline MIC_CALLABLE_MEMBER Complex<long double> mod( const Complex<long double>& x, const Complex<long double>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> log( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> log( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> log( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> floor( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> floor( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> floor( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> ceil( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> ceil( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> ceil( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float arg( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double arg( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double arg( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sin( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sin( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sin( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sinh( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sinh( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sinh( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> cos( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> cos( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> cos( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> cosh( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> cosh( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> cosh( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> tan( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> tan( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> tan( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> atan( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> atan( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> atan( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> copysign( const Complex<float>& x, const Complex<float>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> copysign( const Complex<double>& x, const Complex<double>& y );

    static inline MIC_CALLABLE_MEMBER Complex<long double> copysign( const Complex<long double>& x, const Complex<long double>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float real( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double real( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double real( const Complex<long double>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float imag( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double imag( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double imag( const Complex<long double>& x );

#endif

    using MathReal::sqrt;
    using MathReal::abs;
    using MathReal::conj;
    using MathReal::exp;
    using MathReal::pow;
    using MathReal::mod;
    using MathReal::log;
    using MathReal::floor;
    using MathReal::ceil;
    using MathReal::sin;
    using MathReal::sinh;
    using MathReal::cos;
    using MathReal::cosh;
    using MathReal::tan;
    using MathReal::atan;
    using MathReal::atan2;
    using MathReal::copysign;
    using MathReal::real;
    using MathReal::imag;

};

#ifdef SCAI_COMPLEX_SUPPORTED

// ------------------ Math::sqrt --------------------------------

// Workaround for template function: code for complex abs function as macro
// Reason: template function causes compiler warnings with attribute CUDA_CALLABLE

#define SQRT_FUNCTION_CODE                                                            \
    ValueType x = a.real();                                                           \
    ValueType y = a.imag();                                                           \
                                                                                      \
    ValueType zero = 0;                                                               \
                                                                                      \
    if ( x == zero )                                                                  \
    {                                                                                 \
        ValueType t = MathReal::sqrt( MathReal::abs( y ) / 2 );                       \
        return Complex<ValueType>( t, y < zero ? -t : t );                            \
    }                                                                                 \
    else                                                                              \
    {                                                                                 \
        ValueType t = MathReal::sqrt( 2 * ( Math::abs( a ) + MathReal::abs( x ) ) );  \
        ValueType u = t * ValueType( 0.5 );                                           \
        return x > zero ? Complex<ValueType>( u, y / t ) :                            \
               Complex<ValueType>( MathReal::abs( y ) / t, y < zero ? -u : u );       \
    }    

Complex<float> Math::sqrt( const Complex<float>& a )
{
    typedef float ValueType;

    SQRT_FUNCTION_CODE
}

Complex<double> Math::sqrt( const Complex<double>& a )
{
    typedef double ValueType;

    SQRT_FUNCTION_CODE
}

Complex<long double> Math::sqrt( const Complex<long double>& a )
{   
    typedef long double ValueType;

    SQRT_FUNCTION_CODE
}

// ------------------ Math::abs --------------------------------

// Workaround for template function: code for complex abs function as macro
// Reason: template function causes compiler warnings with attribute CUDA_CALLABLE

#define ABS_FUNCTION_CODE                       \
    ValueType x = a.real();                     \
    ValueType y = a.imag();                     \
    const ValueType ax = MathReal::abs( x );    \
    const ValueType ay = MathReal::abs( y );    \
    const ValueType s  = ax > ay ? ax : ay;     \
                                                \
    if ( s == ValueType( 0 ) )                  \
    {                                           \
        return s;                               \
    }                                           \
                                                \
    x /= s;                                     \
    y /= s;                                     \
                                                \
    return s * MathReal::sqrt( x * x + y * y );

float Math::abs( const Complex<float>& a )
{
    typedef float ValueType;

    ABS_FUNCTION_CODE
}

double Math::abs( const Complex<double>& a )
{
    typedef float ValueType;

    ABS_FUNCTION_CODE
}

long double Math::abs( const Complex<long double>& a )
{
    typedef long double ValueType;

    ABS_FUNCTION_CODE
}

// ------------------ Math::conj --------------------------------

Complex<float> Math::conj( const Complex<float>& a )
{
    return Complex<float>( a.real(), -a.imag() );
}

Complex<double> Math::conj( const Complex<double>& a )
{
    return Complex<double>( a.real(), -a.imag() );
}

Complex<long double> Math::conj( const Complex<long double>& a )
{
    return Complex<long double>( a.real(), -a.imag() );
}

// ------------------ Math::exp --------------------------------

Complex<float> Math::exp( const Complex<float>& a )
{
    float s, c;
    float e = ::expf( a.real() );
    s = ::sinf( a.imag() );
    c = ::cosf( a.imag() );
    return Complex<float>( c * e, s * e );
}
Complex<double> Math::exp( const Complex<double>& a )
{
    double s, c;
    double e = ::exp( a.real() );
    s = ::sin( a.imag() );
    c = ::cos( a.imag() );
    return Complex<double>( c * e, s * e );
}
Complex<long double> Math::exp( const Complex<long double>& a )
{
    long double s, c;
    long double e = ::expl( a.real() );
    s = ::sinl( a.imag() );
    c = ::cosl( a.imag() );
    return Complex<long double>( c * e, s * e );
}

// ------------------ Math::pow --------------------------------

Complex<float> Math::pow( const Complex<float>& base, const Complex<float>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

Complex<double> Math::pow( const Complex<double>& base, const Complex<double>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

Complex<long double> Math::pow( const Complex<long double>& base, const Complex<long double>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

// ------------------ Math::log --------------------------------

Complex<float> Math::log( const Complex<float>& x )
{
    return Complex<float>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

Complex<double> Math::log( const Complex<double>& x )
{
    return Complex<double>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

Complex<long double> Math::log( const Complex<long double>& x )
{
    return Complex<long double>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

// ------------------ Math::floor --------------------------------

Complex<float> Math::floor( const Complex<float>& x )
{
    return Complex<float>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

Complex<double> Math::floor( const Complex<double>& x )
{
    return Complex<double>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

Complex<long double> Math::floor( const Complex<long double>& x )
{
    return Complex<long double>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

// ------------------ Math::mod --------------------------------

Complex<float> Math::mod( const Complex<float>& x, const Complex<float>& y )
{
    return x - floor( x / y ) * y;
}

Complex<double> Math::mod( const Complex<double>& x, const Complex<double>& y )
{
    return x - floor( x / y ) * y;
}

Complex<long double> Math::mod( const Complex<long double>& x, const Complex<long double>& y )
{
    return x - floor( x / y ) * y;
}

// ------------------ Math::ceil --------------------------------

Complex<float> Math::ceil( const Complex<float>& x )
{
    return Complex<float>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

Complex<double> Math::ceil( const Complex<double>& x )
{
    return Complex<double>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

Complex<long double> Math::ceil( const Complex<long double>& x )
{
    return Complex<long double>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

// ------------------ Math::sin --------------------------------

Complex<float> Math::sin( const Complex<float>& x )
{
    return Complex<float>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

Complex<double> Math::sin( const Complex<double>& x )
{
    return Complex<double>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

Complex<long double> Math::sin( const Complex<long double>& x )
{
    return Complex<long double>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

// ------------------ Math::sinh --------------------------------
Complex<float> Math::sinh( const Complex<float>& x )
{   
    return Complex<float>( Math::sinh( x.real() ) * Math::cos( x.imag() ), Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

Complex<double> Math::sinh( const Complex<double>& x )
{   
    return Complex<double>( Math::sinh( x.real() ) * Math::cos( x.imag() ), Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

Complex<long double> Math::sinh( const Complex<long double>& x )
{
    return Complex<long double>( Math::sinh( x.real() ) * Math::cos( x.imag() ), Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

// ------------------ Math::cos --------------------------------
Complex<float> Math::cos( const Complex<float>& x )
{
    return Complex<float>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

Complex<double> Math::cos( const Complex<double>& x )
{
    return Complex<double>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

Complex<long double> Math::cos( const Complex<long double>& x )
{
    return Complex<long double>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

// ------------------ Math::cosh --------------------------------

Complex<float> Math::cosh( const Complex<float>& x )
{
    return Complex<float>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

Complex<double> Math::cosh( const Complex<double>& x )
{
    return Complex<double>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

Complex<long double> Math::cosh( const Complex<long double>& x )
{
    return Complex<long double>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

// ------------------ Math::tan --------------------------------

Complex<float> Math::tan( const Complex<float>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

Complex<double> Math::tan( const Complex<double>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

Complex<long double> Math::tan( const Complex<long double>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

// ------------------ Math::atan --------------------------------

Complex<float> Math::atan( const Complex<float>& x )
{   
    const float r2 = x.real() * x.real();
    const float r1 = float( 1.0 ) - r2 - x.imag() * x.imag();
    
    float num = x.imag() + float( 1.0 );
    float den = x.imag() - float( 1.0 );
    
    num = r2 + num * num;
    den = r2 + den * den;
    
    return Complex<float>( float( 0.5 )  * atan2( float( 2.0 ) * x.real(), r1 ),
                           float( 0.25 ) * log( num / den ) );
}

Complex<double> Math::atan( const Complex<double>& x )
{
    const double r2 = x.real() * x.real();
    const double r1 = double( 1.0 ) - r2 - x.imag() * x.imag();
    
    double num = x.imag() + double( 1.0 );
    double den = x.imag() - double( 1.0 );
    
    num = r2 + num * num;
    den = r2 + den * den;
    
    return Complex<double>( double( 0.5 )  * atan2( double( 2.0 ) * x.real(), r1 ),
                           double( 0.25 ) * log( num / den ) );
}

Complex<long double> Math::atan( const Complex<long double>& x )
{
    const long double r2 = x.real() * x.real();
    const long double r1 = static_cast<long double>( 1.0 ) - r2 - x.imag() * x.imag();

    long double num = x.imag() + static_cast<long double>( 1.0 );
    long double den = x.imag() - static_cast<long double>( 1.0 );

    num = r2 + num * num;
    den = r2 + den * den;

    return Complex<long double>( static_cast<long double>( 0.5 )  * atan2( static_cast<long double>( 2.0 ) * x.real(), r1 ),
                                 static_cast<long double>( 0.25 ) * log( num / den ) );
}

// ------------------ Math::copysign --------------------------------

Complex<float> Math::copysign( const Complex<float>& x, const Complex<float>& y )
{
    return Complex<float>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

Complex<double> Math::copysign( const Complex<double>& x, const Complex<double>& y )
{
    return Complex<double>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

Complex<long double> Math::copysign( const Complex<long double>& x, const Complex<long double>& y )
{
    return Complex<long double>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

// ------------------ Math::arg --------------------------------

float Math::arg( const Complex<float>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

double Math::arg( const Complex<double>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

long double Math::arg( const Complex<long double>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

// ------------------ Math::real --------------------------------

float Math::real( const Complex<float>& a )
{
    return a.real();
}

double Math::real( const Complex<double>& a )
{
    return a.real();
}

long double Math::real( const Complex<long double>& a )
{
    return a.real();
}

// ------------------ Math::imag --------------------------------

float Math::imag( const Complex<float>& a )
{
    return a.imag();
}

double Math::imag( const Complex<double>& a )
{
    return a.imag();
}

long double Math::imag( const Complex<long double>& a )
{
    return a.imag();
}

// ------------------ Math::random ------------------------------

template<>
inline Complex<float> MathReal::random( const unsigned bound )
{
    return Complex<float>( random<float>( bound ), random<float>( bound ) );
}

template<>
inline Complex<double> MathReal::random( const unsigned bound )
{
    return Complex<double>( random<double>( bound ), random<double>( bound ) );
}

template<>
inline Complex<long double> MathReal::random( const unsigned bound )
{
    return Complex<long double>( random<long double>( bound ), random<long double>( bound ) );
}

#endif

} /* end namespace common */

} /* end namespace scai */
