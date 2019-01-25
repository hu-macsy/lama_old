/**
 * @file Complex.hpp
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

#define COMPLEX_CONSTRUCTOR_CUDA( type, method )                                    \
    CUDA_CALLABLE_MEMBER                                                            \
    inline Complex( type value )                                                    \
    {                                                                               \
        method                                                                      \
    }

#define COMPLEX_CONSTRUCTOR_NONCUDA( type, method )                                 \
    inline Complex( type value )                                                    \
    {                                                                               \
        method                                                                      \
    }

#define COMPLEX_CONSTRUCTOR_2_CUDA( type1, type2, method )                          \
    CUDA_CALLABLE_MEMBER                                                            \
    inline Complex( type1 value1, type2 value2 )                                    \
    {                                                                               \
        method                                                                      \
    }

#define COMPLEX_CONSTRUCTOR_2_NONCUDA( type1, type2, method )                       \
    inline Complex( type1 value1, type2 value2 )                                    \
    {                                                                               \
        method                                                                      \
    }

/*
 * For member calculation operators =, +=, -=, *= and /=
 */

#define COMPLEX_OPERATOR_CUDA( op, type, method )                                   \
    CUDA_CALLABLE_MEMBER                                                            \
    inline Complex<ValueType>& op( const type t )                                   \
    {                                                                               \
        method                                                                      \
    }                                                                               \
    CUDA_CALLABLE_MEMBER                                                            \
    inline volatile Complex<ValueType>& op( const type t ) volatile                 \
    {                                                                               \
        method                                                                      \
    }

#define COMPLEX_OPERATOR_NONCUDA( op, type, method )                                \
    inline Complex<ValueType>& op( const type t )                                   \
    {                                                                               \
        method                                                                      \
    }                                                                               \
    inline volatile Complex<ValueType>& op( const type t ) volatile                 \
    {                                                                               \
        method                                                                      \
    }

/*
 * For cast operators
 */

#define COMPLEX_OPERATOR_CAST_CUDA( type, method )                                  \
    CUDA_CALLABLE_MEMBER                                                            \
    operator type() const                                                           \
    {                                                                               \
        method                                                                      \
    }

#define COMPLEX_OPERATOR_CAST_NONCUDA( type, method )                               \
    operator type() const                                                           \
    {                                                                               \
        method                                                                      \
    }

/*
 * For equality operators: ==, !=
 */

#define COMPLEX_OPERATOR_EQUALITY_CUDA(op, sign, connection, type)                            \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline bool op( const Complex<ValueType>& a, const Complex<type>& b )                     \
    {                                                                                         \
        return a.real() sign b.real() connection a.imag() sign b.imag();                      \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline bool op( const Complex<ValueType>& a, const type& b )                              \
    {                                                                                         \
        return a.real() sign b connection a.imag() sign 0;                                    \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline bool op( const type& a, const Complex<ValueType>& b )                              \
    {                                                                                         \
        return a sign b.real() connection 0 sign b.imag();                                    \
    }

#define COMPLEX_OPERATOR_EQUALITY_NONCUDA(op, sign, connection, type)                           \
    template<typename ValueType>                                                                \
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const Complex<type>& b )  \
    {                                                                                           \
        return a.real() sign b.real() connection a.imag() sign b.imag();                        \
    }                                                                                           \
    template<typename ValueType>                                                                \
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const type& b )           \
    {                                                                                           \
        return a.real() sign b connection a.imag() sign 0;                                      \
    }                                                                                           \
    template<typename ValueType>                                                                \
    CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<ValueType>& b )           \
    {                                                                                           \
        return a sign b.real() connection 0 sign b.imag();                                      \
    }

/*
 * For non-member calculation operators: +,-,*,/
 */

#define COMPLEX_OPERATOR_NONMEMBER_CUDA(op, sign, type)                                       \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const Complex<type>& b )       \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<ValueType> op( volatile Complex<ValueType>& a, const Complex<type>& b )    \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )              \
    {                                                                                         \
        Complex<type> x = a;                                                                  \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const type& b )                \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign static_cast<ValueType>(b);                                                     \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<ValueType> op( const type& a, const Complex<ValueType>& b )                \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    CUDA_CALLABLE_MEMBER                                                                      \
    inline Complex<ValueType> op( Complex<ValueType>& a, Complex<type>& b )                   \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }
     
#define COMPLEX_OPERATOR_NONMEMBER_NONCUDA(op, sign, type)                                    \
    template<typename ValueType>                                                              \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const Complex<type>& b )       \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    inline Complex<ValueType> op( volatile Complex<ValueType>& a, const Complex<type>& b )    \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )              \
    {                                                                                         \
        Complex<type> x = a;                                                                  \
        x sign b;                                                                             \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    inline Complex<ValueType> op( const Complex<ValueType>& a, const type& b )                \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign static_cast<ValueType>(b);                                                     \
        return x;                                                                             \
    }                                                                                         \
    template<typename ValueType>                                                              \
    inline Complex<ValueType> op( const type& a, const Complex<ValueType>& b )                \
    {                                                                                         \
        Complex<ValueType> x = a;                                                             \
        x sign b;                                                                             \
        return x;                                                                             \
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

    CUDA_CALLABLE_MEMBER
    explicit Complex();

    template<typename OtherValueType>
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
     * @brief Returns imaginary part of the complex number
     *
     * @return     the imaginary part of a complex number
     */
    CUDA_CALLABLE_MEMBER
    inline ValueType imag( void )
    {
        return mParts[1];
    }

    CUDA_CALLABLE_MEMBER
    inline ValueType imag( void ) const
    {
        return mParts[1];
    }

    CUDA_CALLABLE_MEMBER
    inline ValueType imag( void ) volatile
    {
        return mParts[1];
    }

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
    CUDA_CALLABLE_MEMBER
    inline void imag( const ValueType a )
    {
        mParts[1] = a;
    }

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
    CUDA_CALLABLE_MEMBER
    inline ValueType real( void )
    {
        return mParts[0];
    }

    CUDA_CALLABLE_MEMBER
    inline ValueType real( void ) const
    {
        return mParts[0];
    }

    CUDA_CALLABLE_MEMBER
    inline ValueType real( void ) volatile
    {
        return mParts[0];
    }

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
    CUDA_CALLABLE_MEMBER
    inline void real( const ValueType a )
    {
        mParts[0] = a;
    }

    CUDA_CALLABLE_MEMBER
    inline void real( const ValueType a ) volatile
    {
        mParts[0] = a;
    }

private:
    ValueType mParts[2]; // mParts[0] <-- Real, //mParts[1] <-- Imag
};

template<typename ValueType>
CUDA_CALLABLE_MEMBER
Complex<ValueType>::Complex()
{
    //Constructor have to be empty. Otherwise cuda will output a lot of warnings
}

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

} /* end namespace common */

// define more convenient names

typedef common::Complex<float> ComplexFloat;
typedef common::Complex<double> ComplexDouble;
typedef common::Complex<long double> ComplexLongDouble;

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
#undef COMPLEX_OPERATOR_EQUALITY_CUDA
#undef COMPLEX_OPERATOR_EQUALITY_NONCUDA
#undef COMPLEX_OPERATOR_NONMEMBER_CUDA
#undef COMPLEX_OPERATOR_NONMEMBER_NONCUDA

