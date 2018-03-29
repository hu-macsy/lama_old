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

#include <scai/common/TypeTraits.hpp>
#include <scai/common/cuda/CUDACallable.hpp>

#ifdef SCAI_COMPLEX_SUPPORTED
#include <scai/common/Complex.hpp>
#endif

namespace scai
{

namespace common
{

/** All functions appear whereever possible as template routines.
 *
 *  This is not possible for abs, real, imag 
 */
struct Math 
{
    /** Absolute value function for ValueType
     *
     *  In contrary to the routine of cmath it will be possible to
     *  use always the same name for the routine.
     */

    template<typename ValueType>
    static inline typename TypeTraits<ValueType>::RealType abs( const ValueType& x );

    /*
     * Getter for the real part
     */
    template<typename ValueType>
    static inline typename TypeTraits<ValueType>::RealType real( const ValueType& x );

    /*
     * Getter for the imag part
     */
    template<typename ValueType>
    static inline typename TypeTraits<ValueType>::RealType imag( const ValueType& x );

    /** Square root function for ValueType
     *
     *  In contrary to the routine of cmath it will be possible to
     *  use always the same name for the routine.
     *
     *  \code
     *    ValueType x = sqrt ( y );            // might not work always correctly
     *    ValueType x = Math::sqrt ( y );      // this is guaranteed to work for all arithmetic types
     *  \endcode
     */
    template<typename ValueType>
    static inline ValueType sqrt( const ValueType& x );

    /*
     * Computes the conjugated value of a given value
     */
    template<typename ValueType>
    static inline ValueType conj( const ValueType& x );

    /*
     * Computes the exponential function of a given value
     */
    template<typename ValueType>
    static inline ValueType exp( const ValueType& x );

    /*
     * pow-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType pow( const ValueType& base, const ValueType& exponent );

    /*
     * mod-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType mod( const ValueType& x, const ValueType& y );

    /*
     * log-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType log( const ValueType& x );

    /*
     * floor-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType floor( const ValueType& x );

    /*
     * ceil-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType ceil( const ValueType& x );

    /*
     * sin-function for ValueType
     */

    template<typename ValueType>
    static inline ValueType sin( const ValueType& x );

    /*
     * sinh-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType sinh( const ValueType& x );

    /*
     * cos-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType cos( const ValueType& x );

    /*
     * cosh-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType cosh( const ValueType& x );

    /*
     * tan-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType tan( const ValueType& x );

    /*
     * atan-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType atan( const ValueType& x );

    /*
     * arg function for ValueType, only used for complex
     */
    template<typename ValueType>
    static inline typename TypeTraits<ValueType>::RealType arg( const ValueType& x );

    /*
     * atan2-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType atan2( const ValueType& y, const ValueType& x );

    /*
     * copysign-function for ValueType
     */
    template<typename ValueType>
    static inline ValueType copysign( const ValueType& y, const ValueType& x );

    /*
     * min operation
     */
    template<typename ValueType>
    static inline CUDA_CALLABLE_MEMBER ValueType min( const ValueType& x, const ValueType& y );

    /*
     * max operation
     */
    template<typename ValueType>
    static inline CUDA_CALLABLE_MEMBER ValueType max( const ValueType& x, const ValueType& y );

    /*
     * @brief random value betweem 0 and bound inclusive
     */
    template<typename ValueType>
    static inline ValueType random( unsigned bound );

    /** 
     * @brief initialize seed of random number generator
     */
    static inline void srandom( unsigned int seed );

    /** @brief generate a boolean value by random with a certain ratio.
     *
     *  @param[in] trueRatio specifies probability for true, 0 returns always false, 1 returns always tru
     *  @returns   a boolean value
     */
    static inline bool randomBool( const float trueRatio );
 
    /** @brief return the exponent p for the smallest power that satisfies $2^p \ge n$ */

    template<typename ValueType>
    static inline CUDA_CALLABLE_MEMBER ValueType nextpow2( ValueType n );
};

// -------------------------------- sqrt ----------------------------

template<>
inline CUDA_CALLABLE_MEMBER float Math::sqrt( const float& x )
{
    return ::sqrtf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER double Math::sqrt( const double& x )
{
    return ::sqrt( x );
}

template<>
inline
long double Math::sqrt( const long double& x )
{
    return ::sqrtl( x );
}

// -------------------------------- abs -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
int Math::abs( const int& x )
{
    return ::abs( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
unsigned int Math::abs( const unsigned int& x )
{
    return x;
}

template<>
inline CUDA_CALLABLE_MEMBER
long Math::abs( const long& x )
{
    return ::labs( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
unsigned long Math::abs( const unsigned long& x )
{
    return x;
}

template<>
inline CUDA_CALLABLE_MEMBER
float Math::abs( const float& x )
{
    return ::fabsf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::abs( const double& x )
{
    return ::fabs( x );
}

template<>
inline
long double Math::abs( const long double& x )
{
    return ::fabsl( x );
}

// -------------------------------- conj -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::conj( const float& x )
{
    return x;
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::conj( const double& x )
{
    return x;
}

template<>
inline 
long double
Math::conj( const long double& x )
{
    return x;
}

// -------------------------------- exp -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::exp( const float& x )
{
    return ::expf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::exp( const double& x )
{
    return ::exp( x );
}

template<>
inline 
long double Math::exp( const long double& x )
{
    return ::expl( x );
}

// -------------------------------- pow -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
int Math::pow( const int& base, const int& exponent )
{
    int r = 1;

    for ( int i = 0; i < exponent; ++i )
    {
        r *= base;
    }

    return r;
}

template<>
inline CUDA_CALLABLE_MEMBER
float Math::pow( const float& base, const float& exponent )
{
    return powf( base, exponent );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::pow( const double& base, const double& exponent )
{
    return ::pow( base, exponent );
}

template<>
inline
long double Math::pow( const long double& base, const long double& exponent )
{
    return powl( base, exponent );
}

// -------------------------------- mod -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::mod( const float& x, const float& y )
{
    return fmodf( x, y );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::mod( const double& x, const double& y )
{
    return fmod( x, y );
}

template<>
inline
long double Math::mod( const long double& x, const long double& y )
{
    return fmodl( x, y );
}

// -------------------------------- log -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::log( const float& x )
{
    return logf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::log( const double& x )
{
    return ::log( x );
}

template<>
inline
long double Math::log( const long double& x )
{
    return logl( x );
}

// -------------------------------- floor -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::floor( const float& x )
{
    return floorf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::floor( const double& x )
{
    return ::floor( x );
}

template<>
inline
long double Math::floor( const long double& x )
{
    return floorl( x );
}

// -------------------------------- ceil -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::ceil( const float& x )
{
    return ceilf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::ceil( const double& x )
{
    return ::ceil( x );
}

template<>
inline
long double Math::ceil( const long double& x )
{
    return ceill( x );
}

// -------------------------------- sin -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER 
float Math::sin( const float& x )
{
    return sinf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER 
double Math::sin( const double& x )
{
    return ::sin( x );
}

template<>
inline 
long double Math::sin( const long double& x )
{
    return sinl( x );
}

// -------------------------------- sinh -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER float Math::sinh( const float& x )
{
    return sinhf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER double Math::sinh( const double& x )
{
    return ::sinh( x );
}

template<>
inline long double Math::sinh( const long double& x )
{
    return sinhl( x );
}

// -------------------------------- cos -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::cos( const float& x )
{
    return cosf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::cos( const double& x )
{
    return ::cos( x );
}

template<>
inline 
long double Math::cos( const long double& x )
{
    return cosl( x );
}

// -------------------------------- cosh -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::cosh( const float& x )
{
    return coshf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::cosh( const double& x )
{
    return ::cosh( x );
}

template<>
inline
long double Math::cosh( const long double& x )
{
    return coshl( x );
}

// -------------------------------- tan -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::tan( const float& x )
{
    return tanf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::tan( const double& x )
{
    return ::tan( x );
}

template<>
inline
long double Math::tan( const long double& x )
{
    return tanl( x );
}

// -------------------------------- atan -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::atan( const float& x )
{
    return atanf( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::atan( const double& x )
{
    return ::atan( x );
}

template<>
inline
long double Math::atan( const long double& x )
{
    return atanl( x );
}

// -------------------------------- atan2 -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::atan2( const float& y, const float& x )
{
    return atan2f( y, x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::atan2( const double& y, const double& x )
{
    return ::atan2( y, x );
}

template<>
inline
long double Math::atan2( const long double& y, const long double& x )
{
    return atan2l( y, x );
}

// -------------------------------- copysign -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::copysign( const float& y, const float& x )
{
    return copysignf( y, x );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::copysign( const double& y, const double& x )
{
    return ::copysign( y, x );
}

template<>
inline
long double Math::copysign( const long double& y, const long double& x )
{
    return copysignl( y, x );
}

// -------------------------------- real -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::real( const float& x )
{
    return x;
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::real( const double& x )
{
    return x;
}

template<>
inline
long double Math::real( const long double& x )
{
    return x;
}

// -------------------------------- imag -----------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::imag( const float& )
{
    return 0;
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::imag( const double& )
{
    return 0;
}

template<>
inline
long double Math::imag( const long double& )
{
    return 0;
}

// -------------------------------- min ------------------------------

template<typename ValueType>
inline ValueType Math::min( const ValueType& x, const ValueType& y )
{
    if ( y < x )
    {
        return y;
    }
    else
    {
        return x;
    }
}

// -------------------------------- max ------------------------------

template<typename ValueType>
inline ValueType Math::max( const ValueType& x, const ValueType& y )
{
    if ( x < y )
    {
        return y;
    }
    else
    {
        return x;
    }
}

// -------------------------------- random ---------------------------

template<typename ValueType>
inline ValueType Math::random( const unsigned bound )
{
    return rand() % ( bound + 1 );
}

// Template specializations, must also have the inline attribute

template<>
inline float Math::random( const unsigned bound )
{
    return static_cast<float>( rand() ) / static_cast<float>( RAND_MAX ) * bound;
}

template<>
inline double Math::random( const unsigned bound )
{
    return static_cast<double>( rand() ) / static_cast<double>( RAND_MAX ) * bound;
}

template<>
inline long double Math::random( const unsigned bound )
{
    return static_cast<long double>( rand() ) / static_cast<long double>( RAND_MAX ) * bound;
}

inline void Math::srandom( unsigned int seed )
{
    srand( seed );
}

bool Math::randomBool( const float trueRatio )
{
    if ( trueRatio <= 0.0f )
    {
         return false;
    }
    else if ( trueRatio >= 1.0f )
    {
         return true;
    }
    else
    {
         float x = random<float>( 1 );   // value between 0 and 1, uniformly distributed
         return x < trueRatio;
    }
}

// -------------------------------- nextpow2  ------------------------

template<typename ValueType>
ValueType Math::nextpow2( const ValueType n )
{
    ValueType a = Math::abs( n );
    ValueType p2 = 1;
    ValueType p = 0;

    while ( p2 < a )
    {
        p2 *= 2;
        p  += 1;
    }

    return p;
}

#ifdef SCAI_COMPLEX_SUPPORTED

// ------------------ Math::abs --------------------------------

// Workaround for template function: code for complex abs function as macro
// Reason: template function causes compiler warnings with attribute CUDA_CALLABLE
//
// Note: Make sure to give the same result as real numbers when the
// complex number is purely real or purely imaginary.

#define ABS_FUNCTION_CODE                       \
    ValueType x = a.real();                     \
    ValueType y = a.imag();                     \
    const ValueType ax = Math::abs( x );        \
    const ValueType ay = Math::abs( y );        \
    const ValueType s  = ax > ay ? ax : ay;     \
                                                \
    if ( s == ValueType( 0 ) )                  \
    {                                           \
        return s;                               \
    }                                           \
    else if ( ay == ValueType ( 0 ) )           \
    {                                           \
        return ax;                              \
    }                                           \
    else if ( ax == ValueType ( 0 ) )           \
    {                                           \
        return ay;                              \
    }                                           \
                                                \
    x /= s;                                     \
    y /= s;                                     \
                                                \
    return s * Math::sqrt( x * x + y * y );

template<>
inline CUDA_CALLABLE_MEMBER
float Math::abs( const Complex<float>& a )
{
    typedef float ValueType;

    ABS_FUNCTION_CODE
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::abs( const Complex<double>& a )
{
    typedef double ValueType;

    ABS_FUNCTION_CODE
}

template<>
inline
long double Math::abs( const Complex<long double>& a )
{
    typedef long double ValueType;

    ABS_FUNCTION_CODE
}

// ------------------ Math::sqrt --------------------------------

// Workaround for template function: code for complex function as macro
// Reason: template function causes compiler warnings with attribute CUDA_CALLABLE

#define SQRT_FUNCTION_CODE                                                      \
                                                                                \
    ValueType x = a.real();                                                     \
    ValueType y = a.imag();                                                     \
                                                                                \
    ValueType zero = 0;                                                         \
                                                                                \
    if ( x == zero )                                                            \
    {                                                                           \
        ValueType t = Math::sqrt( Math::abs( y ) / 2 );                         \
        return Complex<ValueType>( t, y < zero ? -t : t );                      \
    }                                                                           \
    else                                                                        \
    {                                                                           \
        ValueType t = Math::sqrt( 2 * ( Math::abs( a ) + Math::abs( x ) ) );    \
        ValueType u = t * ValueType( 0.5 );                                     \
        return x > zero ? Complex<ValueType>( u, y / t ) :                      \
               Complex<ValueType>( Math::abs( y ) / t, y < zero ? -u : u );     \
    }    

template<>
inline CUDA_CALLABLE_MEMBER Complex<float> Math::sqrt( const Complex<float>& a )
{
    typedef float ValueType;

    SQRT_FUNCTION_CODE
}

template<>
inline CUDA_CALLABLE_MEMBER Complex<double> Math::sqrt( const Complex<double>& a )
{
    typedef double ValueType;

    SQRT_FUNCTION_CODE
}

template<>
inline
Complex<long double> Math::sqrt( const Complex<long double>& a )
{   
    typedef long double ValueType;

    SQRT_FUNCTION_CODE
}

// ------------------ Math::conj --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::conj( const Complex<float>& a )
{
    return Complex<float>( a.real(), -a.imag() );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::conj( const Complex<double>& a )
{
    return Complex<double>( a.real(), -a.imag() );
}

template<>
inline // not on CUDA
Complex<long double> Math::conj( const Complex<long double>& a )
{
    return Complex<long double>( a.real(), -a.imag() );
}

// ------------------ Math::exp --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::exp( const Complex<float>& a )
{
    float s, c;
    float e = ::expf( a.real() );
    s = ::sinf( a.imag() );
    c = ::cosf( a.imag() );
    return Complex<float>( c * e, s * e );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::exp( const Complex<double>& a )
{
    double s, c;
    double e = ::exp( a.real() );
    s = ::sin( a.imag() );
    c = ::cos( a.imag() );
    return Complex<double>( c * e, s * e );
}

template<>
inline Complex<long double> Math::exp( const Complex<long double>& a )
{
    long double s, c;
    long double e = ::expl( a.real() );
    s = ::sinl( a.imag() );
    c = ::cosl( a.imag() );
    return Complex<long double>( c * e, s * e );
}

// ------------------ Math::arg ( < log ) ----------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::arg( const Complex<float>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::arg( const Complex<double>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

template<>
inline
long double Math::arg( const Complex<long double>& x )
{
    return Math::atan2( x.imag(), x.real() );
}

// ------------------ Math::log  ( < pow )  --------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::log( const Complex<float>& x )
{
    return Complex<float>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::log( const Complex<double>& x )
{
    return Complex<double>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

template<>
inline
Complex<long double> Math::log( const Complex<long double>& x )
{
    return Complex<long double>( Math::log( Math::abs( x ) ), Math::arg( x ) );
}

// ------------------ Math::pow --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER 
Complex<float> Math::pow( const Complex<float>& base, const Complex<float>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

template<>
inline CUDA_CALLABLE_MEMBER 
Complex<double> Math::pow( const Complex<double>& base, const Complex<double>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

template<>
inline Complex<long double> Math::pow( const Complex<long double>& base, const Complex<long double>& exponent )
{
    return Math::exp( exponent * Math::log( base ) );
}

// ------------------ Math::floor --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::floor( const Complex<float>& x )
{
    return Complex<float>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::floor( const Complex<double>& x )
{
    return Complex<double>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

template<>
inline 
Complex<long double> Math::floor( const Complex<long double>& x )
{
    return Complex<long double>( Math::floor( x.real() ), Math::floor( x.imag() ) );
}

// ------------------ Math::mod --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::mod( const Complex<float>& x, const Complex<float>& y )
{
    return x - floor( x / y ) * y;
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::mod( const Complex<double>& x, const Complex<double>& y )
{
    return x - floor( x / y ) * y;
}

template<>
inline
Complex<long double> Math::mod( const Complex<long double>& x, const Complex<long double>& y )
{
    return x - floor( x / y ) * y;
}

// ------------------ Math::ceil --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::ceil( const Complex<float>& x )
{
    return Complex<float>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::ceil( const Complex<double>& x )
{
    return Complex<double>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

template<>
inline
Complex<long double> Math::ceil( const Complex<long double>& x )
{
    return Complex<long double>( Math::ceil( x.real() ), Math::ceil( x.imag() ) );
}

// ------------------ Math::sin --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::sin( const Complex<float>& x )
{
    return Complex<float>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::sin( const Complex<double>& x )
{
    return Complex<double>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

template<>
inline 
Complex<long double> Math::sin( const Complex<long double>& x )
{
    return Complex<long double>( Math::sin( x.real() ) * Math::cosh( x.imag() ), Math::cos( x.real() ) * Math::sinh( x.imag() ) );
}

// ------------------ Math::sinh --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER 
Complex<float> Math::sinh( const Complex<float>& x )
{   
    return Complex<float>( Math::sinh( x.real() ) * Math::cos( x.imag() ), Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::sinh( const Complex<double>& x )
{   
    return Complex<double>( Math::sinh( x.real() ) * Math::cos( x.imag() ), 
                            Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

template<>
inline
Complex<long double> Math::sinh( const Complex<long double>& x )
{
    return Complex<long double>( Math::sinh( x.real() ) * Math::cos( x.imag() ), 
                                 Math::cosh( x.real() ) * Math::sin( x.imag() ) );
}

// ------------------ Math::cos --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::cos( const Complex<float>& x )
{
    return Complex<float>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::cos( const Complex<double>& x )
{
    return Complex<double>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

template<>
inline 
Complex<long double> Math::cos( const Complex<long double>& x )
{
    return Complex<long double>( Math::cos( x.real() ) * Math::cosh( x.imag() ) , -Math::sin( x.real() ) * Math::sinh( x.imag() ) );
}

// ------------------ Math::cosh --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::cosh( const Complex<float>& x )
{
    return Complex<float>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::cosh( const Complex<double>& x )
{
    return Complex<double>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

template<>
inline 
Complex<long double> Math::cosh( const Complex<long double>& x )
{
    return Complex<long double>( Math::cosh( x.real() ) * Math::cos( x.imag() ) , Math::sinh( x.real() ) * Math::sin( x.imag() ) );
}

// ------------------ Math::tan --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::tan( const Complex<float>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::tan( const Complex<double>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

template<>
inline 
Complex<long double> Math::tan( const Complex<long double>& x )
{
    return Math::sin( x ) / Math::cos( x );
}

// ------------------ Math::atan --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
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

template<>
inline CUDA_CALLABLE_MEMBER
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

template<>
inline
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

template<>
inline CUDA_CALLABLE_MEMBER
Complex<float> Math::copysign( const Complex<float>& x, const Complex<float>& y )
{
    return Complex<float>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

template<>
inline CUDA_CALLABLE_MEMBER
Complex<double> Math::copysign( const Complex<double>& x, const Complex<double>& y )
{
    return Complex<double>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

template<>
inline
Complex<long double> Math::copysign( const Complex<long double>& x, const Complex<long double>& y )
{
    return Complex<long double>( copysign( x.real(), y.real() ), copysign( x.imag(), y.imag() ) );
}

// ------------------ Math::real --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::real( const Complex<float>& a )
{
    return a.real();
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::real( const Complex<double>& a )
{
    return a.real();
}

template<>
inline
long double Math::real( const Complex<long double>& a )
{
    return a.real();
}

// ------------------ Math::imag --------------------------------

template<>
inline CUDA_CALLABLE_MEMBER
float Math::imag( const Complex<float>& a )
{
    return a.imag();
}

template<>
inline CUDA_CALLABLE_MEMBER
double Math::imag( const Complex<double>& a )
{
    return a.imag();
}

template<>
inline
long double Math::imag( const Complex<long double>& a )
{
    return a.imag();
}

// ------------------ Math::random ------------------------------

template<>
inline Complex<float> Math::random( const unsigned bound )
{
    return Complex<float>( random<float>( bound ), random<float>( bound ) );
}

template<>
inline Complex<double> Math::random( const unsigned bound )
{
    return Complex<double>( random<double>( bound ), random<double>( bound ) );
}

template<>
inline Complex<long double> Math::random( const unsigned bound )
{
    return Complex<long double>( random<long double>( bound ), random<long double>( bound ) );
}

#endif

} /* end namespace common */

} /* end namespace scai */
