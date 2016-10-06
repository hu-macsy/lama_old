/**
 * @file Math.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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

#ifdef SCAI_COMPLEX_SUPPORTED
template<typename ValueType> class Complex;
#endif

/** Structure (instead of namespace) that contains the required mathematical functions
 *  for each supported arithmetic type.
 */

struct Math
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

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sqrt( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sqrt( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sqrt( const Complex<long double>& x );
#endif

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

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double abs( const Complex<long double>& x );
#endif

    /*
     * Computes the conjugated value of a given value
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float conj( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double conj( const double& x );

    static inline MIC_CALLABLE_MEMBER long double conj( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> conj( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> conj( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> conj( const Complex<long double>& x );
#endif

    /*
     * Computes the exponential function of a given value
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float exp( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double exp( const double& x );

    static inline MIC_CALLABLE_MEMBER long double exp( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> exp( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> exp( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> exp( const Complex<long double>& x );
#endif

    /*
     * pow-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int pow( const int& base, const int& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float pow( const float& base, const float& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double pow( const double& base, const double& exponent );

    static inline MIC_CALLABLE_MEMBER long double pow( const long double& base, const long double& exponent );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> pow( const Complex<float>& base, const Complex<float>& exponent );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> pow( const Complex<double>& base, const Complex<double>& exponent );

    static inline MIC_CALLABLE_MEMBER Complex<long double> pow( const Complex<long double>& base, const Complex<long double>& exponent );
#endif

    /*
     * log-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float log( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double log( const double& x );

    static inline MIC_CALLABLE_MEMBER long double log( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> log( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> log( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> log( const Complex<long double>& x );
#endif

    /*
     * arg-function for ValueType
     */
#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float arg( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double arg( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double arg( const Complex<long double>& x );
#endif

    /*
     * sin-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sin( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sin( const double& x );

    static inline MIC_CALLABLE_MEMBER long double sin( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sin( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sin( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sin( const Complex<long double>& x );
#endif

    /*
     * sinh-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sinh( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sinh( const double& x );

    static inline MIC_CALLABLE_MEMBER long double sinh( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sinh( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sinh( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> sinh( const Complex<long double>& x );
#endif

    /*
     * cos-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float cos( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double cos( const double& x );

    static inline MIC_CALLABLE_MEMBER long double cos( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> cos( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> cos( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> cos( const Complex<long double>& x );
#endif

    /*
     * cos-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float cosh( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double cosh( const double& x );

    static inline MIC_CALLABLE_MEMBER long double cosh( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> cosh( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> cosh( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> cosh( const Complex<long double>& x );
#endif

    /*
     * tan-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float tan( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double tan( const double& x );

    static inline MIC_CALLABLE_MEMBER long double tan( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> tan( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> tan( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> tan( const Complex<long double>& x );
#endif

    /*
     * atan-function for ValueType
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float atan( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double atan( const double& x );

    static inline MIC_CALLABLE_MEMBER long double atan( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> atan( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> atan( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> atan( const Complex<long double>& x );
#endif

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

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> copysign( const Complex<float>& x, const Complex<float>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> copysign( const Complex<double>& x, const Complex<double>& y );

    static inline MIC_CALLABLE_MEMBER Complex<long double> copysign( const Complex<long double>& x, const Complex<long double>& y );
#endif

    /*
     * Getter for the real part
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float real( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double real( const double& x );

    static inline MIC_CALLABLE_MEMBER long double real( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float real( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double real( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double real( const Complex<long double>& x );
#endif

    /*
     * Getter for the imag part
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float imag( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double imag( const double& x );

    static inline MIC_CALLABLE_MEMBER long double imag( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float imag( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double imag( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER long double imag( const Complex<long double>& x );
#endif

    /*
     * min operation
     */

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int min( const int& x, const int& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned int min( const unsigned int& x, const unsigned int& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long min( const long& x, const long& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned long min( const unsigned long& x, const unsigned long& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float min( const float& x, const float& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double min( const double& x, const double& y );

    static inline MIC_CALLABLE_MEMBER long double min( const long double& x, const long double& y );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> min( const Complex<float>& x, const Complex<float>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> min( const Complex<double>& x, const Complex<double>& y );

    static inline MIC_CALLABLE_MEMBER Complex<long double> min( const Complex<long double>& x, const Complex<long double>& y );
#endif

    /*
     * max operation
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int max( const int& x, const int& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned int max( const unsigned int& x, const unsigned int& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long max( const long& x, const long& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER unsigned long max( const unsigned long& x, const unsigned long& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float max( const float& x, const float& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double max( const double& x, const double& y );

    static inline MIC_CALLABLE_MEMBER long double max( const long double& x, const long double& y );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> max( const Complex<float>& x, const Complex<float>& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> max( const Complex<double>& x, const Complex<double>& y );

    static inline MIC_CALLABLE_MEMBER Complex<long double> max( const Complex<long double>& x, const Complex<long double>& y );
#endif

    /*
     * random value creator
     */
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( int& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( unsigned int& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( unsigned long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( double& x );

    static inline MIC_CALLABLE_MEMBER  void random( long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER  void random( Complex<long double>& x );
#endif
};

// -------------------------------- sqrt ----------------------------

float Math::sqrt( const float& x )
{
    return ::sqrtf( x );
}

double Math::sqrt( const double& x )
{
    return ::sqrt( x );
}

long double Math::sqrt( const long double& x )
{
    return ::sqrtl( x );
}

// -------------------------------- abs -----------------------------

int Math::abs( const int& x )
{
    return ::abs( x );
}

unsigned int Math::abs( const unsigned int& x )
{
    return x;
}

long Math::abs( const long& x )
{
    return ::labs( x );
}

unsigned long Math::abs( const unsigned long& x )
{
    return x;
}

float Math::abs( const float& x )
{
    return ::fabsf( x );
}

double Math::abs( const double& x )
{
    return ::fabs( x );
}

long double Math::abs( const long double& x )
{
    return ::fabsl( x );
}

// -------------------------------- conj -----------------------------

float Math::conj( const float& x )
{
    return x;
}

double Math::conj( const double& x )
{
    return x;
}

long double Math::conj( const long double& x )
{
    return x;
}

// -------------------------------- exp -----------------------------

float Math::exp( const float& x )
{
    return ::expf(x);
}

double Math::exp( const double& x )
{
    return ::exp(x);
}

long double Math::exp( const long double& x )
{
    return ::expl(x);
}

// -------------------------------- pow -----------------------------

int Math::pow( const int& base, const int& exponent )
{
    int r = 1;

    for( int i = 0; i < exponent; ++i )
    {
        r *= base;
    }

    return r;
}

float Math::pow( const float& base, const float& exponent )
{
    return powf( base, exponent );
}

double Math::pow( const double& base, const double& exponent )
{
    return ::pow( base, exponent );
}

long double Math::pow( const long double& base, const long double& exponent )
{
    return powl( base, exponent );
}

// -------------------------------- log -----------------------------

float Math::log( const float& x )
{
    return logf( x );
}

double Math::log( const double& x )
{
    return ::log( x );
}

long double Math::log( const long double& x )
{
    return logl( x );
}

// -------------------------------- sin -----------------------------

float Math::sin( const float& x )
{
    return sinf( x );
}

double Math::sin( const double& x )
{
    return ::sin( x );
}

long double Math::sin( const long double& x )
{
    return sinl( x );
}

// -------------------------------- sinh -----------------------------

float Math::sinh( const float& x )
{
    return sinhf( x );
}

double Math::sinh( const double& x )
{
    return ::sinh( x );
}

long double Math::sinh( const long double& x )
{
    return sinhl( x );
}

// -------------------------------- cos -----------------------------

float Math::cos( const float& x )
{
    return cosf( x );
}

double Math::cos( const double& x )
{
    return ::cos( x );
}

long double Math::cos( const long double& x )
{
    return cosl( x );
}

// -------------------------------- cosh -----------------------------

float Math::cosh( const float& x )
{
    return coshf( x );
}

double Math::cosh( const double& x )
{
    return ::cosh( x );
}

long double Math::cosh( const long double& x )
{
    return coshl( x );
}

// -------------------------------- tan -----------------------------

float Math::tan( const float& x )
{
    return tanf( x );
}

double Math::tan( const double& x )
{
    return ::tan( x );
}

long double Math::tan( const long double& x )
{
    return tanl( x );
}

// -------------------------------- atan -----------------------------

float Math::atan( const float& x )
{
    return atanf( x );
}

double Math::atan( const double& x )
{
    return ::atan( x );
}

long double Math::atan( const long double& x )
{
    return atanl( x );
}

// -------------------------------- atan2 -----------------------------

float Math::atan2( const float& y, const float& x )
{
    return atan2f( y, x );
}

double Math::atan2( const double& y, const double& x )
{
    return ::atan2( y, x );
}

long double Math::atan2( const long double& y, const long double& x )
{
    return atan2l( y, x );
}

// -------------------------------- copysign -----------------------------

float Math::copysign( const float& y, const float& x )
{
    return copysignf( y, x );
}

double Math::copysign( const double& y, const double& x )
{
    return ::copysign( y, x );
}

long double Math::copysign( const long double& y, const long double& x )
{
    return copysignl( y, x );
}

// -------------------------------- real -----------------------------

float Math::real( const float& x )
{
    return x;
}

double Math::real( const double& x )
{
    return x;
}

long double Math::real( const long double& x )
{
    return x;
}

// -------------------------------- imag -----------------------------

float Math::imag( const float& )
{
    return 0;
}

double Math::imag( const double& )
{
    return 0;
}

long double Math::imag( const long double& )
{
    return 0;
}

// -------------------------------- min ------------------------------

int Math::min( const int& x, const int& y )
{
    return y < x ? y : x;
}

unsigned int Math::min( const unsigned int& x, const unsigned int& y )
{
    return y < x ? y : x;
}

long Math::min( const long& x, const long& y )
{
    return y < x ? y : x;
}

unsigned long Math::min( const unsigned long& x, const unsigned long& y )
{
    return y < x ? y : x;
}

float Math::min( const float& x, const float& y )
{
    return y < x ? y : x;
}

double Math::min( const double& x, const double& y )
{
    return y < x ? y : x;
}

long double Math::min( const long double& x, const long double& y )
{
    return y < x ? y : x;
}

// -------------------------------- max ------------------------------

int Math::max( const int& x, const int& y )
{
    return x > y ? x : y;
}

unsigned int Math::max( const unsigned int& x, const unsigned int& y )
{
    return x > y ? x : y;
}

long Math::max( const long& x, const long& y )
{
    return x > y ? x : y;
}

unsigned long Math::max( const unsigned long& x, const unsigned long& y )
{
    return x > y ? x : y;
}

float Math::max( const float& x, const float& y )
{
    return x > y ? x : y;
}

double Math::max( const double& x, const double& y )
{
    return x > y ? x : y;
}

long double Math::max( const long double& x, const long double& y )
{
    return x > y ? x : y;
}

// -------------------------------- random ---------------------------

void Math::random( int& x )
{
    x = rand();
}

void Math::random( unsigned int& x )
{
    x = rand();
}

void Math::random( long& x )
{
    x = rand();
}

void Math::random( unsigned long& x )
{
    x = rand();
}

void Math::random( float& x )
{
    x = 1 - static_cast<float>( rand() ) / static_cast<float>( RAND_MAX / 2 );
}

void Math::random( double& x )
{
    x = 1 - static_cast<double>( rand() ) / static_cast<double>( RAND_MAX / 2 );
}

void Math::random( long double& x )
{
    x = 1 - static_cast<long double>( rand() ) / static_cast<long double>( RAND_MAX / 2 );
}

} /* end namespace common */

} /* end namespace scai */
