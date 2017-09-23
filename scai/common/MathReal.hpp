/**
 * @file MathReal.hpp
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

#include <scai/common/mic/MICCallable.hpp>
#include <scai/common/cuda/CUDACallable.hpp>

#include <cmath>
#include <cstdlib>

namespace scai
{

namespace common
{

/** Structure (instead of namespace) that contains the required mathematical functions
 *  for each supported arithmetic type.
 *  
 *  MathReal only for real data types, will be extended later to Math with Complex types
 */
struct MathReal
{
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
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sqrt( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sqrt( const double& x );

    static inline MIC_CALLABLE_MEMBER long double sqrt( const long double& x );

    /** Absolute value function for ValueType
     *
     *  In contrary to the routine of cmath it will be possible to
     *  use always the same name for the routine.
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int abs( const int& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned int abs( const unsigned int& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long abs( const long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned long abs( const unsigned long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const double& x );

    static inline MIC_CALLABLE_MEMBER long double abs( const long double& x );

    /*
     * Computes the conjugated value of a given value
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float conj( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double conj( const double& x );

    static inline MIC_CALLABLE_MEMBER long double conj( const long double& x );

    /*
     * Computes the exponential function of a given value
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float exp( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double exp( const double& x );

    static inline MIC_CALLABLE_MEMBER long double exp( const long double& x );

    /*
     * pow-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int pow( const int& base, const int& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float pow( const float& base, const float& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double pow( const double& base, const double& exponent );

    static inline MIC_CALLABLE_MEMBER long double pow( const long double& base, const long double& exponent );

    /*
     * mod-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float mod( const float& x, const float& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double mod( const double& x, const double& y );

    static inline MIC_CALLABLE_MEMBER long double mod( const long double& x, const long double& y );

    /*
     * log-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float log( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double log( const double& x );

    static inline MIC_CALLABLE_MEMBER long double log( const long double& x );

    /*
     * floor-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float floor( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double floor( const double& x );

    static inline MIC_CALLABLE_MEMBER long double floor( const long double& x );

    /*
     * ceil-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float ceil( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double ceil( const double& x );

    static inline MIC_CALLABLE_MEMBER long double ceil( const long double& x );

    /*
     * sin-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sin( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sin( const double& x );

    static inline MIC_CALLABLE_MEMBER long double sin( const long double& x );

    /*
     * sinh-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sinh( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sinh( const double& x );

    static inline MIC_CALLABLE_MEMBER long double sinh( const long double& x );

    /*
     * cos-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float cos( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double cos( const double& x );

    static inline MIC_CALLABLE_MEMBER long double cos( const long double& x );


    /*
     * cos-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float cosh( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double cosh( const double& x );

    static inline MIC_CALLABLE_MEMBER long double cosh( const long double& x );

    /*
     * tan-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float tan( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double tan( const double& x );

    static inline MIC_CALLABLE_MEMBER long double tan( const long double& x );

    /*
     * atan-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float atan( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double atan( const double& x );

    static inline MIC_CALLABLE_MEMBER long double atan( const long double& x );

    /*
     * atan2-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float atan2( const float& y, const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double atan2( const double& y, const double& x );

    static inline MIC_CALLABLE_MEMBER long double atan2( const long double& y, const long double& x );

    /*
     * copysign-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float copysign( const float& x, const float& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double copysign( const double& x, const double& y );

    static inline MIC_CALLABLE_MEMBER long double copysign( const long double& x, const long double& y );

    /*
     * Getter for the real part
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float real( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double real( const double& x );

    static inline MIC_CALLABLE_MEMBER long double real( const long double& x );

    /*
     * Getter for the imag part
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float imag( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double imag( const double& x );

    static inline MIC_CALLABLE_MEMBER long double imag( const long double& x );

    /*
     * min operation
     */
    template<typename ValueType>
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER ValueType min( const ValueType& x, const ValueType& y );

    /*
     * max operation
     */
    template<typename ValueType>
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER ValueType max( const ValueType& x, const ValueType& y );

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
};

// -------------------------------- sqrt ----------------------------

float MathReal::sqrt( const float& x )
{
    return ::sqrtf( x );
}

double MathReal::sqrt( const double& x )
{
    return ::sqrt( x );
}

long double MathReal::sqrt( const long double& x )
{
    return ::sqrtl( x );
}

// -------------------------------- abs -----------------------------

int MathReal::abs( const int& x )
{
    return ::abs( x );
}

unsigned int MathReal::abs( const unsigned int& x )
{
    return x;
}

long MathReal::abs( const long& x )
{
    return ::labs( x );
}

unsigned long MathReal::abs( const unsigned long& x )
{
    return x;
}

float MathReal::abs( const float& x )
{
    return ::fabsf( x );
}

double MathReal::abs( const double& x )
{
    return ::fabs( x );
}

long double MathReal::abs( const long double& x )
{
    return ::fabsl( x );
}

// -------------------------------- conj -----------------------------

float MathReal::conj( const float& x )
{
    return x;
}

double MathReal::conj( const double& x )
{
    return x;
}

long double MathReal::conj( const long double& x )
{
    return x;
}

// -------------------------------- exp -----------------------------

float MathReal::exp( const float& x )
{
    return ::expf( x );
}

double MathReal::exp( const double& x )
{
    return ::exp( x );
}

long double MathReal::exp( const long double& x )
{
    return ::expl( x );
}

// -------------------------------- pow -----------------------------

int MathReal::pow( const int& base, const int& exponent )
{
    int r = 1;

    for ( int i = 0; i < exponent; ++i )
    {
        r *= base;
    }

    return r;
}

float MathReal::pow( const float& base, const float& exponent )
{
    return powf( base, exponent );
}

double MathReal::pow( const double& base, const double& exponent )
{
    return ::pow( base, exponent );
}

long double MathReal::pow( const long double& base, const long double& exponent )
{
    return powl( base, exponent );
}

// -------------------------------- mod -----------------------------

float MathReal::mod( const float& x, const float& y )
{
    return fmodf( x, y );
}

double MathReal::mod( const double& x, const double& y )
{
    return fmod( x, y );
}

long double MathReal::mod( const long double& x, const long double& y )
{
    return fmodl( x, y );
}

// -------------------------------- log -----------------------------

float MathReal::log( const float& x )
{
    return logf( x );
}

double MathReal::log( const double& x )
{
    return ::log( x );
}

long double MathReal::log( const long double& x )
{
    return logl( x );
}

// -------------------------------- floor -----------------------------

float MathReal::floor( const float& x )
{
    return floorf( x );
}

double MathReal::floor( const double& x )
{
    return ::floor( x );
}

long double MathReal::floor( const long double& x )
{
    return floorl( x );
}

// -------------------------------- ceil -----------------------------

float MathReal::ceil( const float& x )
{
    return ceilf( x );
}

double MathReal::ceil( const double& x )
{
    return ::ceil( x );
}

long double MathReal::ceil( const long double& x )
{
    return ceill( x );
}

// -------------------------------- sin -----------------------------

float MathReal::sin( const float& x )
{
    return sinf( x );
}

double MathReal::sin( const double& x )
{
    return ::sin( x );
}

long double MathReal::sin( const long double& x )
{
    return sinl( x );
}

// -------------------------------- sinh -----------------------------

float MathReal::sinh( const float& x )
{
    return sinhf( x );
}

double MathReal::sinh( const double& x )
{
    return ::sinh( x );
}

long double MathReal::sinh( const long double& x )
{
    return sinhl( x );
}

// -------------------------------- cos -----------------------------

float MathReal::cos( const float& x )
{
    return cosf( x );
}

double MathReal::cos( const double& x )
{
    return ::cos( x );
}

long double MathReal::cos( const long double& x )
{
    return cosl( x );
}

// -------------------------------- cosh -----------------------------

float MathReal::cosh( const float& x )
{
    return coshf( x );
}

double MathReal::cosh( const double& x )
{
    return ::cosh( x );
}

long double MathReal::cosh( const long double& x )
{
    return coshl( x );
}

// -------------------------------- tan -----------------------------

float MathReal::tan( const float& x )
{
    return tanf( x );
}

double MathReal::tan( const double& x )
{
    return ::tan( x );
}

long double MathReal::tan( const long double& x )
{
    return tanl( x );
}

// -------------------------------- atan -----------------------------

float MathReal::atan( const float& x )
{
    return atanf( x );
}

double MathReal::atan( const double& x )
{
    return ::atan( x );
}

long double MathReal::atan( const long double& x )
{
    return atanl( x );
}

// -------------------------------- atan2 -----------------------------

float MathReal::atan2( const float& y, const float& x )
{
    return atan2f( y, x );
}

double MathReal::atan2( const double& y, const double& x )
{
    return ::atan2( y, x );
}

long double MathReal::atan2( const long double& y, const long double& x )
{
    return atan2l( y, x );
}

// -------------------------------- copysign -----------------------------

float MathReal::copysign( const float& y, const float& x )
{
    return copysignf( y, x );
}

double MathReal::copysign( const double& y, const double& x )
{
    return ::copysign( y, x );
}

long double MathReal::copysign( const long double& y, const long double& x )
{
    return copysignl( y, x );
}

// -------------------------------- real -----------------------------

float MathReal::real( const float& x )
{
    return x;
}

double MathReal::real( const double& x )
{
    return x;
}

long double MathReal::real( const long double& x )
{
    return x;
}

// -------------------------------- imag -----------------------------

float MathReal::imag( const float& )
{
    return 0;
}

double MathReal::imag( const double& )
{
    return 0;
}

long double MathReal::imag( const long double& )
{
    return 0;
}

// -------------------------------- min ------------------------------

template<typename ValueType>
inline ValueType MathReal::min( const ValueType& x, const ValueType& y )
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
inline ValueType MathReal::max( const ValueType& x, const ValueType& y )
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
inline ValueType MathReal::random( const unsigned bound )
{
    return rand() % ( bound + 1 );
}

// Template specializations, must also have the inline attribute

template<>
inline float MathReal::random( const unsigned bound )
{
    return static_cast<float>( rand() ) / static_cast<float>( RAND_MAX ) * bound;
}

template<>
inline double MathReal::random( const unsigned bound )
{
    return static_cast<double>( rand() ) / static_cast<double>( RAND_MAX ) * bound;
}

template<>
inline long double MathReal::random( const unsigned bound )
{
    return static_cast<long double>( rand() ) / static_cast<long double>( RAND_MAX ) * bound;
}

inline void MathReal::srandom( unsigned int seed )
{
    srand( seed );
}

bool MathReal::randomBool( const float trueRatio )
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

} /* end namespace common */

} /* end namespace scai */
