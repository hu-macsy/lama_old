/**
 * @file Complex.hpp
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
 * @brief Definition of template class Complex
 * @author Eric Schricker
 * @date 22.01.2014
 */

#pragma once

/** Enable preprocessor flag to tell other code that ComplexZZZ is supported */
#define SCAI_COMPLEX_SUPPORTED

// local library
#include <scai/common/config.hpp>
#include <scai/common/cuda/CUDACallable.hpp>
#include <scai/common/mic/MICCallable.hpp>
#include <scai/common/Math.hpp>

// std
#include <sstream>
#include <algorithm>

#ifdef __CUDACC__
#include <cuComplex.h>
//#include <math.h>
#include <thrust/device_reference.h>
#else
//#include <cmath>
#endif

/*
 * Macros used for building the operators, providing the functionality
 * Every macro is defined with the suffix REAL and COMPLEX. Each is
 * providing another functionality.
 */

#define COMPLEX_SET_REAL                                                                                \
    real( static_cast<ValueType>( t ) );                                                                \
    imag( static_cast<ValueType>( 0 ) );                                                                \
    return *this;

#define COMPLEX_SET_COMPLEX                                                                             \
    real( static_cast<ValueType>( t.real() ) );                                                         \
    imag( static_cast<ValueType>( t.imag() ) );                                                         \
    return *this;

#define COMPLEX_OP_REAL( op )                                                                           \
    real( real() op static_cast<ValueType>( t ) );                                                      \
    return *this;

#define COMPLEX_OP_COMPLEX( op )                                                                        \
    real( real() op static_cast<ValueType>( t.real() ) );                                               \
    imag( imag() op static_cast<ValueType>( t.imag() ) );                                               \
    return *this;

#define COMPLEX_MULTIPLY_REAL                                                                           \
    real( real() * static_cast<ValueType>( t ) );                                                       \
    imag( imag() * static_cast<ValueType>( t ) );                                                       \
    return *this;

#define COMPLEX_MULTIPLY_COMPLEX                                                                        \
    ValueType a = real();                                                                               \
    ValueType b = imag();                                                                               \
    ValueType c = static_cast<ValueType>( t.real() );                                                   \
    ValueType d = static_cast<ValueType>( t.imag() );                                                   \
    real( a * c - b * d );                                                                              \
    imag( a * d + b * c );                                                                              \
    return *this;

#define COMPLEX_DIVIDE_REAL                                                                             \
    real( real() / static_cast<ValueType>( t ) );                                                       \
    imag( imag() / static_cast<ValueType>( t ) );                                                       \
    return *this;

#define COMPLEX_DIVIDE_COMPLEX                                                                          \
    ValueType a = real();                                                                               \
    ValueType b = imag();                                                                               \
    ValueType c = static_cast<ValueType>( t.real() );                                                   \
    ValueType d = static_cast<ValueType>( t.imag() );                                                   \
    ValueType c2d2 = ( c * c + d * d);                                                                  \
    real( ( a * c + b * d ) / ( c2d2 ) );                                                               \
    imag( ( b * c - a * d ) / ( c2d2 ) );                                                               \
    return *this;

#define COMPLEX_CAST_REAL(type)                                                                         \
    return static_cast<type>( real() );

#define COMPLEX_CAST_COMPLEX(type)                                                                      \
    return Complex<type>( static_cast<type>( real() ), static_cast<type>( imag() ) );

#define COMPLEX_CONSTRUCTOR_REAL                                                                        \
    this->real( static_cast<ValueType>( value ) );                                                      \
    this->imag( static_cast<ValueType>( 0 ) );

#define COMPLEX_CONSTRUCTOR_REAL2                                                                       \
    this->real( static_cast<ValueType>( value1 ) );                                                     \
    this->imag( static_cast<ValueType>( value2 ) );

#define COMPLEX_CONSTRUCTOR_COMPLEX                                                                     \
    this->real( static_cast<ValueType>( value.real() ) );                                               \
    this->imag( static_cast<ValueType>( value.imag() ) );

/*
 * Macros which build the operators, the macros above will be used for the functionality
 * The macros with the suffix CUDA are getting in an cuda environment the prefix __host__ __device__
 * The macros with the suffix NONCUDA are getting no prefix
 *
 * For creating constructors
 *
 */

#define COMPLEX_CONSTRUCTOR_CUDA( type, method )                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline Complex( type value )                                                                        \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_CONSTRUCTOR_NONCUDA( type, method )                                                     \
    inline Complex( type value )                                                                        \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_CONSTRUCTOR_2_CUDA( type1, type2, method )                                              \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline Complex( type1 value1, type2 value2 )                                                        \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_CONSTRUCTOR_2_NONCUDA( type1, type2, method )                                           \
    inline Complex( type1 value1, type2 value2 )                                                        \
    {                                                                                                   \
        method                                                                                          \
    }

/*
 * For member calculation operators =, +=, -=, *= and /=
 */

#define COMPLEX_OPERATOR_CUDA( op, type, method )                                                       \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline Complex<ValueType>& op( const type t )                                                       \
    {                                                                                                   \
        method                                                                                          \
    }                                                                                                   \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline volatile Complex<ValueType>& op( const type t ) volatile                                     \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_OPERATOR_NONCUDA( op, type, method )                                                    \
    inline Complex<ValueType>& op( const type t )                                                       \
    {                                                                                                   \
        method                                                                                          \
    }                                                                                                   \
    inline volatile Complex<ValueType>& op( const type t ) volatile                                     \
    {                                                                                                   \
        method                                                                                          \
    }

/*
 * For cast operators
 */

#define COMPLEX_OPERATOR_CAST_CUDA( type, method )                                                      \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    operator type() const                                                                               \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_OPERATOR_CAST_NONCUDA( type, method )                                                   \
    operator type() const                                                                               \
    {                                                                                                   \
        method                                                                                          \
    }

/*
 * For comparison operators: <, >, <=, >=
 *
 */

#define COMPLEX_OPERATOR_COMPARISON_CUDA(op, sign, type)                                                \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const Complex<ValueType>& a, const Complex<type>& b )                               \
    {                                                                                                   \
        return a.metrikCuda() sign static_cast<ValueType>(b.metrikCuda());                              \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const Complex<ValueType>& a, const type& b )                                        \
    {                                                                                                   \
        return a.metrikCuda() sign static_cast<ValueType>(b);                                           \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const type& a, const Complex<ValueType>& b )                                        \
    {                                                                                                   \
        return a sign static_cast<type>(b.metrikCuda());                                                \
    }

#define COMPLEX_OPERATOR_COMPARISON_NONCUDA(op, sign, type)                                             \
    template<typename ValueType>                                                                        \
    inline bool op( const Complex<ValueType>& a, const Complex<type>& b )                               \
    {                                                                                                   \
        return a.metrikHost() sign static_cast<ValueType>(b.metrikHost());                              \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    inline bool op( const Complex<ValueType>& a, const type& b )                                        \
    {                                                                                                   \
        return a.metrikHost() sign static_cast<ValueType>(b);                                           \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    inline bool op( const type& a, const Complex<ValueType>& b )                                        \
    {                                                                                                   \
        return a sign static_cast<type>(b.metrikHost());                                                \
    }

/*
 * For equality operators: ==, !=
 */

#define COMPLEX_OPERATOR_EQUALITY_CUDA(op, sign, connection, type)                                      \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const Complex<ValueType>& a, const Complex<type>& b )                               \
    {                                                                                                   \
        return a.real() sign b.real() connection a.imag() sign b.imag();                                \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const Complex<ValueType>& a, const type& b )                                        \
    {                                                                                                   \
        return a.real() sign b connection a.imag() sign 0;                                              \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    MIC_CALLABLE_MEMBER                                                                                 \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline bool op( const type& a, const Complex<ValueType>& b )                                        \
    {                                                                                                   \
        return a sign b.real() connection 0 sign b.imag();                                              \
    }

#define COMPLEX_OPERATOR_EQUALITY_NONCUDA(op, sign, connection, type)                                   \
    template<typename ValueType>                                                                        \
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const Complex<type>& b )          \
    {                                                                                                   \
        return a.real() sign b.real() connection a.imag() sign b.imag();                                \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const type& b )                   \
    {                                                                                                   \
        return a.real() sign b connection a.imag() sign 0;                                              \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<ValueType>& b )                   \
    {                                                                                                   \
        return a sign b.real() connection 0 sign b.imag();                                              \
    }

/*
 * For non-member calculation operators: +,-,*,/
 */

#define COMPLEX_OPERATOR_NONMEMBER_CUDA(op, sign, type)                                                         \
    template<typename ValueType>                                                                                \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const Complex<type>& b )                         \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<ValueType> op( volatile Complex<ValueType>& a, const Complex<type>& b )                      \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )                                \
    {                                                                                                           \
        Complex<type> x = a;                                                                                    \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const type& b )                                  \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign static_cast<ValueType>(b);                                                                       \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<ValueType> op( const type& a, const Complex<ValueType>& b )                                  \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    MIC_CALLABLE_MEMBER                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                                        \
    inline Complex<ValueType> op( Complex<ValueType>& a, Complex<type>& b )                                     \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \

#define COMPLEX_OPERATOR_NONMEMBER_NONCUDA(op, sign, type)                                              \
    template<typename ValueType>                                                                        \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const Complex<type>& b )                 \
    {                                                                                                   \
        Complex<ValueType> x = a;                                                                       \
        x sign b;                                                                                       \
        return x;                                                                                       \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    inline Complex<ValueType> op( volatile Complex<ValueType>& a, const Complex<type>& b )              \
    {                                                                                                   \
        Complex<ValueType> x = a;                                                                       \
        x sign b;                                                                                       \
        return x;                                                                                       \
    }                                                                                                   \
    inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )                        \
    {                                                                                                   \
        Complex<type> x = a;                                                                            \
        x sign b;                                                                                       \
        return x;                                                                                       \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const type& b )                          \
    {                                                                                                   \
        Complex<ValueType> x = a;                                                                       \
        x sign static_cast<ValueType>(b);                                                               \
        return x;                                                                                       \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    inline Complex<ValueType> op( const type& a, const Complex<ValueType>& b )                          \
    {                                                                                                   \
        Complex<ValueType> x = a;                                                                       \
        x sign b;                                                                                       \
        return x;                                                                                       \
    }

    namespace scai
    {

    namespace common
    {


    /**
     * @brief The class Complex represents complex numbers.
     */

    template<typename ValueType>
    class COMMON_DLL_IMPORTEXPORT Complex
    {
    public:
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        explicit Complex();

        template<typename OtherValueType>
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        explicit inline Complex( const OtherValueType value );

        /*
         * @brief Constructs a complex representing the passed real
         *
         * @param[in] real the real part this complex should represent
         */
        COMPLEX_CONSTRUCTOR_CUDA( const char, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const short, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const int, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const unsigned int, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const long, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const unsigned long, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const float, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_CUDA( const double, COMPLEX_CONSTRUCTOR_REAL )
        COMPLEX_CONSTRUCTOR_NONCUDA( const long double, COMPLEX_CONSTRUCTOR_REAL )

        /*
         * @brief Constructs a complex representing the passed real and imaginary part
         *
         * @param[in] real the real part this complex should represent
         * @param[in] imag the imaginary part this complex should represent
         */
        COMPLEX_CONSTRUCTOR_2_CUDA( const int, const int, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const int, const long, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const int, const float, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const int, const double, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const int, const long double, COMPLEX_CONSTRUCTOR_REAL2 )

        COMPLEX_CONSTRUCTOR_2_CUDA( const long, const int, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const long, const long, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const long, const float, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const long, const double, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long, const long double, COMPLEX_CONSTRUCTOR_REAL2 )

        COMPLEX_CONSTRUCTOR_2_CUDA( const float, const int, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const float, const long, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const float, const float, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const float, const double, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const float, const long double, COMPLEX_CONSTRUCTOR_REAL2 )

        COMPLEX_CONSTRUCTOR_2_CUDA( const double, const int, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const double, const long, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const double, const float, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_CUDA( const double, const double, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const double, const long double, COMPLEX_CONSTRUCTOR_REAL2 )

        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long double, const int, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long double, const long, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long double, const float, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long double, const double, COMPLEX_CONSTRUCTOR_REAL2 )
        COMPLEX_CONSTRUCTOR_2_NONCUDA( const long double, const long double, COMPLEX_CONSTRUCTOR_REAL2 )

#ifdef __CUDACC__
        template<typename OtherValueType>
        CUDA_CALLABLE_MEMBER inline Complex( const thrust::device_reference<Complex<OtherValueType> >& value )
        {
            *this = static_cast<Complex<ValueType> >( value );
        }
#endif

        /*
         * @brief Constructs a scalar representing the passed complex value.
         *
         * @param[in]   value the value this complex should represent
         */
        COMPLEX_CONSTRUCTOR_CUDA( const Complex<int>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( volatile const Complex<int>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( const Complex<long>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( volatile const Complex<long>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( const Complex<float>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( volatile const Complex<float>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( const Complex<double>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_CUDA( volatile const Complex<double>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_NONCUDA( const Complex<long double>&, COMPLEX_CONSTRUCTOR_COMPLEX )
        COMPLEX_CONSTRUCTOR_NONCUDA( volatile const Complex<long double>&, COMPLEX_CONSTRUCTOR_COMPLEX )

        /*
         * Creation of the member operators through macros
         */

        /*
         * @brief Binary operator
         */
        COMPLEX_OPERATOR_CUDA( operator=, int, COMPLEX_SET_REAL )
        COMPLEX_OPERATOR_CUDA( operator=, long, COMPLEX_SET_REAL )
        COMPLEX_OPERATOR_CUDA( operator=, float, COMPLEX_SET_REAL )
        COMPLEX_OPERATOR_CUDA( operator=, double, COMPLEX_SET_REAL )
        COMPLEX_OPERATOR_NONCUDA( operator=, long double, COMPLEX_SET_REAL )
        COMPLEX_OPERATOR_CUDA( operator=, Complex<int>, COMPLEX_SET_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator=, Complex<long>, COMPLEX_SET_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator=, Complex<float>, COMPLEX_SET_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator=, Complex<double>, COMPLEX_SET_COMPLEX )
        COMPLEX_OPERATOR_NONCUDA( operator=, Complex<long double>, COMPLEX_SET_COMPLEX )
#ifdef __CUDACC__
        template<typename OtherValueType>
        CUDA_CALLABLE_MEMBER Complex<ValueType>& operator=( const thrust::device_reference<Complex<OtherValueType> >& value )
        {
            Complex<OtherValueType> x = value;
            real( static_cast<ValueType>( x.real() ) );
            imag( static_cast<ValueType>( x.imag() ) );
            return *this;
        }
#endif

        COMPLEX_OPERATOR_CUDA( operator+=, int, COMPLEX_OP_REAL( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, long, COMPLEX_OP_REAL( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, float, COMPLEX_OP_REAL( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, double, COMPLEX_OP_REAL( + ) )
        COMPLEX_OPERATOR_NONCUDA( operator+=, long double, COMPLEX_OP_REAL( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, Complex<int>, COMPLEX_OP_COMPLEX( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, Complex<long>, COMPLEX_OP_COMPLEX( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, Complex<float>, COMPLEX_OP_COMPLEX( + ) )
        COMPLEX_OPERATOR_CUDA( operator+=, Complex<double>, COMPLEX_OP_COMPLEX( + ) )
        COMPLEX_OPERATOR_NONCUDA( operator+=, Complex<long double>, COMPLEX_OP_COMPLEX( + ) )

        COMPLEX_OPERATOR_CUDA( operator-=, int, COMPLEX_OP_REAL( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, long, COMPLEX_OP_REAL( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, float, COMPLEX_OP_REAL( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, double, COMPLEX_OP_REAL( - ) )
        COMPLEX_OPERATOR_NONCUDA( operator-=, long double, COMPLEX_OP_REAL( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, Complex<int>, COMPLEX_OP_COMPLEX( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, Complex<long>, COMPLEX_OP_COMPLEX( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, Complex<float>, COMPLEX_OP_COMPLEX( - ) )
        COMPLEX_OPERATOR_CUDA( operator-=, Complex<double>, COMPLEX_OP_COMPLEX( - ) )
        COMPLEX_OPERATOR_NONCUDA( operator-=, Complex<long double>, COMPLEX_OP_COMPLEX( - ) )

        COMPLEX_OPERATOR_CUDA( operator*=, int, COMPLEX_MULTIPLY_REAL )
        COMPLEX_OPERATOR_CUDA( operator*=, long, COMPLEX_MULTIPLY_REAL )
        COMPLEX_OPERATOR_CUDA( operator*=, float, COMPLEX_MULTIPLY_REAL )
        COMPLEX_OPERATOR_CUDA( operator*=, double, COMPLEX_MULTIPLY_REAL )
        COMPLEX_OPERATOR_NONCUDA( operator*=, long double, COMPLEX_MULTIPLY_REAL )
        COMPLEX_OPERATOR_CUDA( operator*=, Complex<int>, COMPLEX_MULTIPLY_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator*=, Complex<long>, COMPLEX_MULTIPLY_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator*=, Complex<float>, COMPLEX_MULTIPLY_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator*=, Complex<double>, COMPLEX_MULTIPLY_COMPLEX )
        COMPLEX_OPERATOR_NONCUDA( operator*=, Complex<long double>, COMPLEX_MULTIPLY_COMPLEX )

        COMPLEX_OPERATOR_CUDA( operator/=, int, COMPLEX_DIVIDE_REAL )
        COMPLEX_OPERATOR_CUDA( operator/=, long, COMPLEX_DIVIDE_REAL )
        COMPLEX_OPERATOR_CUDA( operator/=, float, COMPLEX_DIVIDE_REAL )
        COMPLEX_OPERATOR_CUDA( operator/=, double, COMPLEX_DIVIDE_REAL )
        COMPLEX_OPERATOR_NONCUDA( operator/=, long double, COMPLEX_DIVIDE_REAL )
        COMPLEX_OPERATOR_CUDA( operator/=, Complex<int>, COMPLEX_DIVIDE_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator/=, Complex<long>, COMPLEX_DIVIDE_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator/=, Complex<float>, COMPLEX_DIVIDE_COMPLEX )
        COMPLEX_OPERATOR_CUDA( operator/=, Complex<double>, COMPLEX_DIVIDE_COMPLEX )
        COMPLEX_OPERATOR_NONCUDA( operator/=, Complex<long double>, COMPLEX_DIVIDE_COMPLEX )

        /*
         * @brief Casts operator
         */
        COMPLEX_OPERATOR_CAST_CUDA( int, COMPLEX_CAST_REAL( int ) )
        COMPLEX_OPERATOR_CAST_CUDA( unsigned int, COMPLEX_CAST_REAL( unsigned int ) )
        COMPLEX_OPERATOR_CAST_CUDA( long, COMPLEX_CAST_REAL( long ) )
        COMPLEX_OPERATOR_CAST_CUDA( unsigned long, COMPLEX_CAST_REAL( unsigned long ) )
        COMPLEX_OPERATOR_CAST_CUDA( float, COMPLEX_CAST_REAL( float ) )
        COMPLEX_OPERATOR_CAST_CUDA( double, COMPLEX_CAST_REAL( double ) )
        COMPLEX_OPERATOR_CAST_NONCUDA( long double, COMPLEX_CAST_REAL( long double ) )

        /*
         * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is callable from CUDA space.
         *
         * @return      the metric of a complex number
         */
        inline MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        ValueType metrikCuda( void ) const;

        /*
         * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is not callable from CUDA space and used for long double operations.
         *
         * @return      the metric of a complex number
         */
        ValueType metrikHost( void ) const;

        /*
         * @brief Returns imaginary part of the complex number
         *
         * @return     the imaginary part of a complex number
         */
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType imag( void )
        {
            return mParts[1];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType imag( void ) const
        {
            return mParts[1];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType imag( void ) volatile
        {
            return mParts[1];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType imag( void ) const volatile
        {
            return mParts[1];
        }

        /*
         * @brief Sets imaginary part of the complex number
         *
         * @param[in] a    the new imaginary part of the complex number
         */
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline void imag( const ValueType a )
        {
            mParts[1] = a;
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline void imag( const ValueType a ) volatile
        {
            mParts[1] = a;
        }

        /**
         * @brief Returns real part of the complex number
         *
         * @return     the real part of a complex number
         */
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType real( void )
        {
            return mParts[0];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType real( void ) const
        {
            return mParts[0];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType real( void ) volatile
        {
            return mParts[0];
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline ValueType real( void ) const volatile
        {
            return mParts[0];
        }

        /**
         * @brief Sets real part of the complex number
         *
         * @param[in] a    the new real part of the complex number
         */
        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline void real( const ValueType a )
        {
            mParts[0] = a;
        }

        MIC_CALLABLE_MEMBER
        CUDA_CALLABLE_MEMBER
        inline void real( const ValueType a ) volatile
        {
            mParts[0] = a;
        }

    private:
        ValueType mParts[2]; // mParts[0] <-- Real, //mParts[1] <-- Imag
    };

    template<typename ValueType>
    MIC_CALLABLE_MEMBER
    CUDA_CALLABLE_MEMBER
    Complex<ValueType>::Complex()
    {
        //Constructor have to be empty. Otherwise cuda will output a lot of warnings
    }

    template<typename ValueType>
    MIC_CALLABLE_MEMBER
    CUDA_CALLABLE_MEMBER
    ValueType Complex<ValueType>::metrikCuda( void ) const
    {
        if( imag() == ValueType( 0 ) )
        {
            // saves time and keeps precision, e.g. for eps0

            return Math::abs( real() );
        }
        else
        {
            return Math::sqrt( real() * real() + imag() * imag() );
        }
    }

    template<typename ValueType>
    ValueType Complex<ValueType>::metrikHost( void ) const
    {
        if ( imag() == ValueType( 0 ) )
        {
            // saves time and keeps precision, e.g. for eps0

            return Math::abs( real() );
        }
        else
        {
            return Math::sqrt( real() * real() + imag() * imag() );
        }
    }

    /*
     * @brief Check if a is lower than b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is lower than b
     */
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, < , int )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, < , long )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, < , float )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, < , double )
    COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator<, < , long double )

    /*
     * @brief Check if a is greater than b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is greater than b
     */
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, > , int )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, > , long )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, > , float )
    COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, > , double )
    COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator>, > , long double )

    /*
     * @brief Check if a is lower than b or equal to b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is lower than b or equal to b
     */
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <= , int )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <= , long )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <= , float )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <= , double )
//COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator<=, <= , long double )

    /*
     * @brief Check if a is greater than b or equal to b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is greater than b or equal to b
     */
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >= , int )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >= , long )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >= , float )
//COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >= , double )
//COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator>=, >= , long double )

    /*
     * @brief Check equality of a and b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is equal to b
     */
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, == , && , int )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, == , && , long )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, == , && , float )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, == , && , double )
    COMPLEX_OPERATOR_EQUALITY_NONCUDA( operator==, == , && , long double )

    /*
     * @brief Check inequality of a and b.
     *
     * @param[in] a     the 1st Complex to compare this to.
     * @param[in] b     the 2nd Complex to compare this to.
     * @return          if a is unequal to b
     */
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, != , || , int )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, != , || , long )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, != , || , float )
    COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, != , || , double )
    COMPLEX_OPERATOR_EQUALITY_NONCUDA( operator!=, != , || , long double )

    template<typename ValueType>
    MIC_CALLABLE_MEMBER
    CUDA_CALLABLE_MEMBER
    inline Complex<ValueType> operator-( const Complex<ValueType>& a )
    {
        Complex<ValueType> x;
        x.real( -a.real() );
        x.imag( -a.imag() );
        return x;
    }

    /*
     * @brief Add Complex a with Complex b
     *
     * @param[in] a     1st Complex.
     * @param[in] b     2nd Complex.
     * @return          sum
     */
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, += , int )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, += , long )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, += , float )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, += , double )
    COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator+, += , long double )

    /*
     * @brief Subtract Complex a with Complex b
     *
     * @param[in] a     1st Complex.
     * @param[in] b     2nd Complex.
     * @return          difference
     */
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -= , int )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -= , long )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -= , float )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -= , double )
    COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator-, -= , long double )

    /*
     * @brief Multiply Complex a with Complex b
     *
     * @param[in] a     1st Complex.
     * @param[in] b     2nd Complex.
     * @return          product
     */
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *= , int )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *= , long )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *= , float )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *= , double )
    COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator*, *= , long double )

    /*
     * @brief Divide Complex a with Complex b
     *
     * @param[in] a     1st Complex.
     * @param[in] b     2nd Complex.
     * @return          quotient
     */
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /= , int )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /= , long )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /= , float )
    COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /= , double )
    COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator/, /= , long double )

    template<typename ValueType, typename InputType1, typename InputType2>
    std::basic_istream<InputType1, InputType2>&
    operator>>( std::basic_istream<InputType1, InputType2>& input, Complex<ValueType>& x )
    {
        ValueType real = 0;
        ValueType imag = 0;
        input >> real;
        if ( !input.fail() )
        {
            input >> imag;
            // a fail on second argument is ignored
            input.clear();
        }
        x.real( real );
        x.imag( imag );
        return input;
    }

    template<typename ValueType>
    std::ostream& operator<<( std::ostream& stream, const Complex<ValueType>& object )
    {
        if ( object.imag() == 0 )
        {
            stream << object.real();
        }
        else
        {
            stream << object.real() << " " << object.imag();
        }

        return stream;
    }

    /*

    template<>
    long double Complex<long double>::metrikCuda( void ) const
    {
        return Math::sqrt( real() * real() + imag() * imag() );
    }

    */

// ------------------ Math::sqrt --------------------------------

    Complex<float> Math::sqrt( const Complex<float>& a )
    {
        float x = a.real();
        float y = a.imag();

        if ( x == 0.0f )
        {
            float t = Math::sqrt( Math::abs( y ) / 2 );
            return Complex<float>( t, y < 0.0f ? -t : t );
        }
        else
        {
            float t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
            float u = t / 2;
            return x > 0.0f ? Complex<float>( u, y / t ) :
                   Complex<float>( Math::abs( y ) / t, y < 0.0f ? -u : u );
        }
    }

    Complex<double> Math::sqrt( const Complex<double>& a )
    {
        double x = a.real();
        double y = a.imag();

        if ( x == 0.0 )
        {
            double t = Math::sqrt( Math::abs( y ) / 2 );
            return Complex<double>( t, y < 0.0 ? -t : t );
        }
        else
        {
            double t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
            double u = t / 2;
            return x > 0.0 ? Complex<double>( u, y / t ) :
                   Complex<double>( Math::abs( y ) / t, y < 0.0 ? -u : u );
        }
    }

    Complex<long double> Math::sqrt( const Complex<long double>& a )
    {
        long double x = a.real();
        long double y = a.imag();

        if ( x == 0.0l )
        {
            long double t = Math::sqrt( Math::abs( y ) / 2 );
            return Complex<long double>( t, y < 0.0l ? -t : t );
        }
        else
        {
            long double t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
            long double u = t / 2;
            return x > 0.0l ? Complex<long double>( u, y / t ) :
                   Complex<long double>( Math::abs( y ) / t, y < 0.0l ? -u : u );
        }
    }

// ------------------ Math::abs --------------------------------
    float Math::abs( const Complex<float>& a )
    {
        float x = a.real();
        float y = a.imag();
        const float ax = Math::abs( x );
        const float ay = Math::abs( y );
        const float s = ax > ay ? ax : ay;

        if ( s == 0.0f )
        {
            return s;
        }

        x /= s;
        y /= s;
        return s * Math::sqrt( x * x + y * y );
    }

    double Math::abs( const Complex<double>& a )
    {
        double x = a.real();
        double y = a.imag();
        const double ax = Math::abs( x );
        const double ay = Math::abs( y );
        const double s = ax > ay ? ax : ay;

        if ( s == 0.0 )
        {
            return s;
        }

        x /= s;
        y /= s;
        return s * Math::sqrt( x * x + y * y );
    }

    long double Math::abs( const Complex<long double>& a )
    {
        long double x = a.real();
        long double y = a.imag();
        const long double ax = Math::abs( x );
        const long double ay = Math::abs( y );
        const long double s = ax > ay ? ax : ay;

        if ( s == 0.0l )
        {
            return s;
        }

        x /= s;
        y /= s;
        return s * Math::sqrt( x * x + y * y );
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
        return Complex<float>( Math::log( Math::abs(x) ), Math::arg(x) );
    }

    Complex<double> Math::log( const Complex<double>& x )
    {
        return Complex<double>( Math::log( Math::abs(x) ), Math::arg(x) );
    }

    Complex<long double> Math::log( const Complex<long double>& x )
    {
        return Complex<long double>( Math::log( Math::abs(x) ), Math::arg(x) );
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
        return Math::sin(x) / Math::cos(x);
    }

    Complex<double> Math::tan( const Complex<double>& x )
    {
        return Math::sin(x) / Math::cos(x);
    }

    Complex<long double> Math::tan( const Complex<long double>& x )
    {
        return Math::sin(x) / Math::cos(x);
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

        return Complex<float>( float( 0.5 )  * atan2( float( 2.0 ) * x.real(), r1),
                               float( 0.25 ) * log( num / den) );
    }

    Complex<double> Math::atan( const Complex<double>& x )
    {
        const double r2 = x.real() * x.real();
        const double r1 = double( 1.0 ) - r2 - x.imag() * x.imag();

        double num = x.imag() + double( 1.0 );
        double den = x.imag() - double( 1.0 );

        num = r2 + num * num;
        den = r2 + den * den;

        return Complex<double>( double( 0.5 )  * atan2( double( 2.0 ) * x.real(), r1),
                               double( 0.25 ) * log( num / den) );
    }

    Complex<long double> Math::atan( const Complex<long double>& x )
    {
        const long double r2 = x.real() * x.real();
        const long double r1 = static_cast<long double>( 1.0 ) - r2 - x.imag() * x.imag();

        long double num = x.imag() + static_cast<long double>( 1.0 );
        long double den = x.imag() - static_cast<long double>( 1.0 );

        num = r2 + num * num;
        den = r2 + den * den;

        return Complex<long double>( static_cast<long double>( 0.5 )  * atan2( static_cast<long double>( 2.0 ) * x.real(), r1),
                               static_cast<long double>( 0.25 ) * log( num / den) );
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

// ------------------ Math::min ---------------------------------
    Complex<float> Math::min( const Complex<float>& x, const Complex<float>& y )
    {
        return Math::abs( x ) < Math::abs( y ) ? x : y;
    }

    Complex<double> Math::min( const Complex<double>& x, const Complex<double>& y )
    {
        return Math::abs( x ) < Math::abs( y ) ? x : y;
    }

    Complex<long double> Math::min( const Complex<long double>& x, const Complex<long double>& y )
    {
        return Math::abs( x ) < Math::abs( y ) ? x : y;
    }

// ------------------ Math::max ---------------------------------
    Complex<float> Math::max( const Complex<float>& x, const Complex<float>& y )
    {
        return Math::abs( x ) > Math::abs( y ) ? x : y;
    }

    Complex<double> Math::max( const Complex<double>& x, const Complex<double>& y )
    {
        return Math::abs( x ) > Math::abs( y ) ? x : y;
    }

    Complex<long double> Math::max( const Complex<long double>& x, const Complex<long double>& y )
    {
        return Math::abs( x ) > Math::abs( y ) ? x : y;
    }

// ------------------ Math::random ------------------------------
    void Math::random( Complex<float>& x )
    {
        float val;
        Math::random( val );
        x.real( val );
        Math::random( val );
        x.imag( val );
    }

    void Math::random( Complex<double>& x )
    {
        double val;
        Math::random( val );
        x.real( val );
        Math::random( val );
        x.imag( val );
    }

    void Math::random( Complex<long double>& x )
    {
        long double val;
        Math::random( val );
        x.real( val );
        Math::random( val );
        x.imag( val );
    }

    } /* end namespace common */

    } /* end namespace scai */

#undef COMPLEX_SET_REAL
#undef COMPLEX_SET_COMPLEX
#undef COMPLEX_OP_REAL
#undef COMPLEX_OP_COMPLEX
#undef COMPLEX_MULTIPLY_REAL
#undef COMPLEX_MULTIPLY_COMPLEX
#undef COMPLEX_DIVIDE_REAL
#undef COMPLEX_DIVIDE_COMPLEX
#undef COMPLEX_CAST_REAL
#undef COMPLEX_CAST_COMPLEX
#undef COMPLEX_CONSTRUCTOR_REAL
#undef COMPLEX_CONSTRUCTOR_REAL2
#undef COMPLEX_CONSTRUCTOR_COMPLEX
#undef COMPLEX_CONSTRUCTOR_CUDA
#undef COMPLEX_CONSTRUCTOR_2_NONCUDA
#undef COMPLEX_CONSTRUCTOR_2_CUDA
#undef COMPLEX_CONSTRUCTOR_NONCUDA
#undef COMPLEX_OPERATOR_CUDA
#undef COMPLEX_OPERATOR_NONCUDA
#undef COMPLEX_OPERATOR_CAST_CUDA
#undef COMPLEX_OPERATOR_CAST_NONCUDA
#undef COMPLEX_OPERATOR_COMPARISON_CUDA
#undef COMPLEX_OPERATOR_COMPARISON_NONCUDA
#undef COMPLEX_OPERATOR_EQUALITY_CUDA
#undef COMPLEX_OPERATOR_EQUALITY_NONCUDA
#undef COMPLEX_OPERATOR_NONMEMBER_CUDA
#undef COMPLEX_OPERATOR_NONMEMBER_NONCUDA

// define more convenient names

    typedef scai::common::Complex<float> ComplexFloat;
    typedef scai::common::Complex<double> ComplexDouble;
    typedef scai::common::Complex<long double> ComplexLongDouble;
