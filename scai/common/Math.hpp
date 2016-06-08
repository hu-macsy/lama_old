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

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long abs( const long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long long abs( const long long& x );

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
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int conj( const int& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long conj( const long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long long conj( const long long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float conj( const float& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double conj( const double& x );

    static inline MIC_CALLABLE_MEMBER long double conj( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> conj( const Complex<float>& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> conj( const Complex<double>& x );

    static inline MIC_CALLABLE_MEMBER Complex<long double> conj( const Complex<long double>& x );
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

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long min( const long& x, const long& y );

    static inline MIC_CALLABLE_MEMBER long long min( const long long& x, const long long& y );

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

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long max( const long& x, const long& y );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long long max( const long long& x, const long long& y );

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

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( long& x );

    static inline MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER void random( long long& x );

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

long Math::abs( const long& x )
{
    return ::labs( x );
}

long long Math::abs( const long long& x )
{
    return ::llabs( x );
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

int Math::conj( const int& x )
{
    return x;
}

long Math::conj( const long& x )
{
    return x;
}

long long Math::conj( const long long& x )
{
    return x;
}

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
    return x < y ? x : y;
}

long Math::min( const long& x, const long& y )
{
    return x < y ? x : y;
}

long long Math::min( const long long& x, const long long& y )
{
    return x < y ? x : y;
}

float Math::min( const float& x, const float& y )
{
    return x < y ? x : y;
}

double Math::min( const double& x, const double& y )
{
    return x < y ? x : y;
}

long double Math::min( const long double& x, const long double& y )
{
    return x < y ? x : y;
}

// -------------------------------- max ------------------------------

int Math::max( const int& x, const int& y )
{
    return x > y ? x : y;
}

long Math::max( const long& x, const long& y )
{
    return x > y ? x : y;
}

long long Math::max( const long long& x, const long long& y )
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

void Math::random( long& x )
{
    x = rand();
}

void Math::random( long long& x )
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
