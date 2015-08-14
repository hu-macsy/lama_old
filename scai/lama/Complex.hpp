/*
 * @file Complex.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Complex.hpp
 * @author Eric Schricker
 * @date 22.01.2014
 * @since 1.1.0
 */

#pragma once

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__ __host__
#include <cuComplex.h>
#include <math.h>
#include <thrust/device_reference.h>
#else
#define CUDA_CALLABLE_MEMBER
#include <cmath>
#endif

#include <scai/lama/exception/LAMAAssert.hpp>

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

#define COMPLEX_PLUS_REAL                                                                               \
    real( real() + static_cast<ValueType>( t ) );                                                       \
    return *this;

#define COMPLEX_PLUS_COMPLEX                                                                            \
    real( real() + static_cast<ValueType>( t.real() ) );                                                \
    imag( imag() + static_cast<ValueType>( t.imag() ) );                                                \
    return *this;

#define COMPLEX_MINUS_REAL                                                                              \
    real( real() - static_cast<ValueType>( t ) );                                                       \
    return *this;

#define COMPLEX_MINUS_COMPLEX                                                                           \
    real( real() - static_cast<ValueType>( t.real() ) );                                                \
    imag( imag() - static_cast<ValueType>( t.imag() ) );                                                \
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
    return static_cast<type>( metrikCuda() );

#define COMPLEX_CAST_COMPLEX(type)                                                                      \
    return Complex<type>( static_cast<type>( real() ), static_cast<type>( imag() ) );

#define COMPLEX_CONSTRUCTOR_REAL                                                                        \
    this->real( static_cast<ValueType>( value ) );                                                      \
    this->imag( static_cast<ValueType>( 0 ) );

/*
 * Macros which build the operators, the macros above will be used for the functionality
 * The macros with the suffix CUDA are getting in an cuda environment the prefix __host__ __device__
 * The macros with the suffix NONCUDA are getting no prefix
 *
 * For creating constructors
 *
 */

#define COMPLEX_CONSTRUCTOR_CUDA( type, method )                                                        \
    CUDA_CALLABLE_MEMBER                                                                                \
    inline Complex( const type value )                                                                  \
    {                                                                                                   \
        method                                                                                          \
    }

#define COMPLEX_CONSTRUCTOR_NONCUDA( type, method )                                                     \
                                                                                                        \
    inline Complex( const type value )                                                                  \
    {                                                                                                   \
        method                                                                                          \
    }

/*
 * For member calculation operators =, +=, -=, *= and /=
 */

#define COMPLEX_OPERATOR_CUDA( op, type, method )                                                       \
    CUDA_CALLABLE_MEMBER inline Complex<ValueType>& op( const type t )                                  \
    {                                                                                                   \
        method                                                                                          \
    }                                                                                                   \
    CUDA_CALLABLE_MEMBER inline volatile Complex<ValueType>& op( const type t ) volatile                \
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
    CUDA_CALLABLE_MEMBER operator type() const                                                          \
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
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const Complex<type>& b )          \
    {                                                                                                   \
        return a.metrikCuda() sign static_cast<ValueType>(b.metrikCuda());                              \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    CUDA_CALLABLE_MEMBER inline bool op( const Complex<ValueType>& a, const type& b )                   \
    {                                                                                                   \
        return a.metrikCuda() sign static_cast<ValueType>(b);                                           \
    }                                                                                                   \
    template<typename ValueType>                                                                        \
    CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<ValueType>& b )                   \
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
    CUDA_CALLABLE_MEMBER inline Complex<ValueType> op( const Complex<ValueType>& a, const Complex<type>& b )    \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    CUDA_CALLABLE_MEMBER inline Complex<ValueType> op( volatile Complex<ValueType>& a, const Complex<type>& b ) \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    CUDA_CALLABLE_MEMBER inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )           \
    {                                                                                                           \
        Complex<type> x = a;                                                                                    \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    CUDA_CALLABLE_MEMBER inline Complex<ValueType> op( const Complex<ValueType>& a, const type& b )             \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign static_cast<ValueType>(b);                                                                       \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    CUDA_CALLABLE_MEMBER inline Complex<ValueType> op( const type& a, const Complex<ValueType>& b )             \
    {                                                                                                           \
        Complex<ValueType> x = a;                                                                               \
        x sign b;                                                                                               \
        return x;                                                                                               \
    }                                                                                                           \
    template<typename ValueType>                                                                                \
    CUDA_CALLABLE_MEMBER inline Complex<ValueType> op( Complex<ValueType>& a, Complex<type>& b )                \
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

/**
 * @brief The class Complex represents complex numbers.
 */
namespace lama
{
using std::sqrt;
using std::abs;

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Complex
{
public:
    CUDA_CALLABLE_MEMBER
    explicit Complex();

    template<typename OtherValueType>
    CUDA_CALLABLE_MEMBER explicit inline Complex( const OtherValueType value );

    /**
     * @brief Constructs a complex representing the passed real
     *
     * @param[in] real the real part this complex should represent
     */
    COMPLEX_CONSTRUCTOR_CUDA(int, COMPLEX_CONSTRUCTOR_REAL)

    /**
     * @brief Constructs a complex representing the passed real
     *
     * @param[in] real the real part this complex should represent
     */
    COMPLEX_CONSTRUCTOR_CUDA(long, COMPLEX_CONSTRUCTOR_REAL)

    /**
     * @brief Constructs a complex representing the passed real
     *
     * @param[in] real the real part this complex should represent
     */
    COMPLEX_CONSTRUCTOR_CUDA(float, COMPLEX_CONSTRUCTOR_REAL)

    /**
     * @brief Constructs a complex representing the passed real
     *
     * @param[in] real the real part this complex should represent
     */
    COMPLEX_CONSTRUCTOR_CUDA(double, COMPLEX_CONSTRUCTOR_REAL)

    /**
     * @brief Constructs a complex representing the passed real
     *
     * @param[in] real the real part this complex should represent
     */
    COMPLEX_CONSTRUCTOR_NONCUDA(long double, COMPLEX_CONSTRUCTOR_REAL)

    /**
     * @brief Constructs a complex representing the passed real and imaginary part
     *
     * @param[in] real the real part this complex should represent
     * @param[in] imag the imaginary part this complex should represent
     */
    template<typename OtherValueType1,typename OtherValueType2>
    CUDA_CALLABLE_MEMBER inline Complex( const OtherValueType1 real, const OtherValueType2 imag )
    {
        this->real( static_cast<ValueType>( real ) );
        this->imag( static_cast<ValueType>( imag ) );
    }

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @param[in]   value the value this complex should represent
     */
    template<typename OtherValueType>
    CUDA_CALLABLE_MEMBER inline Complex( const Complex<OtherValueType>& value )
    {
        this->real( static_cast<ValueType>( value.real() ) );
        this->imag( static_cast<ValueType>( value.imag() ) );
    }

#ifdef __CUDACC__
    template<typename OtherValueType>
    CUDA_CALLABLE_MEMBER inline Complex( const thrust::device_reference<Complex<OtherValueType> >& value )
    {
        *this = static_cast<Complex<ValueType> >(value);
        //this->real( static_cast<ValueType>( (value()).real() ) );
        //this->imag( static_cast<ValueType>( (value()).imag() ) );
    }
#endif

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @param[in]   value the value this complex should represent
     */
    CUDA_CALLABLE_MEMBER
    inline Complex( volatile const Complex<ValueType>& value )
    {
        this->real( value.real() );
        this->imag( value.imag() );
    }

    /*
     * Creation of the member operators through macros
     */

    /**
     * @brief Binary operator
     */
    COMPLEX_OPERATOR_CUDA(operator=, int, COMPLEX_SET_REAL)COMPLEX_OPERATOR_CUDA(operator=, long, COMPLEX_SET_REAL)
    COMPLEX_OPERATOR_CUDA(operator=, float, COMPLEX_SET_REAL)
    COMPLEX_OPERATOR_CUDA(operator=, double, COMPLEX_SET_REAL)
    COMPLEX_OPERATOR_NONCUDA(operator=, long double, COMPLEX_SET_REAL)
    COMPLEX_OPERATOR_CUDA(operator=, Complex<int>, COMPLEX_SET_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator=, Complex<long>, COMPLEX_SET_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator=, Complex<float>, COMPLEX_SET_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator=, Complex<double>, COMPLEX_SET_COMPLEX)
    COMPLEX_OPERATOR_NONCUDA(operator=, Complex<long double>, COMPLEX_SET_COMPLEX)

    COMPLEX_OPERATOR_CUDA(operator+=, int, COMPLEX_PLUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator+=, long, COMPLEX_PLUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator+=, float, COMPLEX_PLUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator+=, double, COMPLEX_PLUS_REAL)
    COMPLEX_OPERATOR_NONCUDA(operator+=, long double, COMPLEX_PLUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator+=, Complex<int>, COMPLEX_PLUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator+=, Complex<long>, COMPLEX_PLUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator+=, Complex<float>, COMPLEX_PLUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator+=, Complex<double>, COMPLEX_PLUS_COMPLEX)
    COMPLEX_OPERATOR_NONCUDA(operator+=, Complex<long double>, COMPLEX_PLUS_COMPLEX)

    COMPLEX_OPERATOR_CUDA(operator-=, int, COMPLEX_MINUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator-=, long, COMPLEX_MINUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator-=, float, COMPLEX_MINUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator-=, double, COMPLEX_MINUS_REAL)
    COMPLEX_OPERATOR_NONCUDA(operator-=, long double, COMPLEX_MINUS_REAL)
    COMPLEX_OPERATOR_CUDA(operator-=, Complex<int>, COMPLEX_MINUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator-=, Complex<long>, COMPLEX_MINUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator-=, Complex<float>, COMPLEX_MINUS_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator-=, Complex<double>, COMPLEX_MINUS_COMPLEX)
    COMPLEX_OPERATOR_NONCUDA(operator-=, Complex<long double>, COMPLEX_MINUS_COMPLEX)

    COMPLEX_OPERATOR_CUDA(operator*=, int, COMPLEX_MULTIPLY_REAL)
    COMPLEX_OPERATOR_CUDA(operator*=, long, COMPLEX_MULTIPLY_REAL)
    COMPLEX_OPERATOR_CUDA(operator*=, float, COMPLEX_MULTIPLY_REAL)
    COMPLEX_OPERATOR_CUDA(operator*=, double, COMPLEX_MULTIPLY_REAL)
    COMPLEX_OPERATOR_NONCUDA(operator*=, long double, COMPLEX_MULTIPLY_REAL)
    COMPLEX_OPERATOR_CUDA(operator*=, Complex<int>, COMPLEX_MULTIPLY_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator*=, Complex<long>, COMPLEX_MULTIPLY_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator*=, Complex<float>, COMPLEX_MULTIPLY_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator*=, Complex<double>, COMPLEX_MULTIPLY_COMPLEX)
    COMPLEX_OPERATOR_NONCUDA(operator*=, Complex<long double>, COMPLEX_MULTIPLY_COMPLEX)

    COMPLEX_OPERATOR_CUDA(operator/=, int, COMPLEX_DIVIDE_REAL)
    COMPLEX_OPERATOR_CUDA(operator/=, long, COMPLEX_DIVIDE_REAL)
    COMPLEX_OPERATOR_CUDA(operator/=, float, COMPLEX_DIVIDE_REAL)
    COMPLEX_OPERATOR_CUDA(operator/=, double, COMPLEX_DIVIDE_REAL)
    COMPLEX_OPERATOR_NONCUDA(operator/=, long double, COMPLEX_DIVIDE_REAL)
    COMPLEX_OPERATOR_CUDA(operator/=, Complex<int>, COMPLEX_DIVIDE_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator/=, Complex<long>, COMPLEX_DIVIDE_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator/=, Complex<float>, COMPLEX_DIVIDE_COMPLEX)
    COMPLEX_OPERATOR_CUDA(operator/=, Complex<double>, COMPLEX_DIVIDE_COMPLEX)
    COMPLEX_OPERATOR_NONCUDA(operator/=, Complex<long double>, COMPLEX_DIVIDE_COMPLEX)

    /**
     * @brief Casts operator
     */
    COMPLEX_OPERATOR_CAST_CUDA(int, COMPLEX_CAST_REAL(int))
    COMPLEX_OPERATOR_CAST_CUDA(long, COMPLEX_CAST_REAL(long))
    COMPLEX_OPERATOR_CAST_CUDA(float, COMPLEX_CAST_REAL(float))
    COMPLEX_OPERATOR_CAST_CUDA(double, COMPLEX_CAST_REAL(double))
    COMPLEX_OPERATOR_CAST_NONCUDA(long double, COMPLEX_CAST_REAL(long double))

    template<typename OtherValueType>
    CUDA_CALLABLE_MEMBER
    operator Complex<OtherValueType>() const
    {
        return Complex<OtherValueType>(static_cast<OtherValueType>(real()), static_cast<OtherValueType>(imag()));
    }

    /**
     * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is callable from CUDA space.
     *
     * @return      the metric of a complex number
     */
    CUDA_CALLABLE_MEMBER ValueType metrikCuda( void ) const;

    /**
     * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is not callable from CUDA space and used for long double operations.
     *
     * @return      the metric of a complex number
     */
    ValueType metrikHost( void ) const;

    /**
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

    /**
     * @brief Sets imaginary part of the complex number
     *
     * @param[in]     the new imaginary part of the complex number
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
     * @param[in]     the new real part of the complex number
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
CUDA_CALLABLE_MEMBER Complex<ValueType>::Complex()
{
    //Constructor have to be empty. Otherwise cuda will output a lot of warnings
}

template<typename ValueType>
CUDA_CALLABLE_MEMBER ValueType Complex<ValueType>::metrikCuda( void ) const
{
    return sqrt( real() * real() + imag() * imag() );
}

template<typename ValueType>
ValueType Complex<ValueType>::metrikHost( void ) const
{
    return sqrt( real() * real() + imag() * imag() );
}

template<typename ValueType>
CUDA_CALLABLE_MEMBER inline ValueType abs( const Complex<ValueType>& a )
{
    return ( sqrt( a.real() * a.real() + a.imag() * a.imag() ) );
}

/*
 * long double version must not be CUDA_CALLABLE_MEMBER
 */
inline long double abs( const Complex<long double>& a )
{
    return ( sqrt( a.real() * a.real() + a.imag() * a.imag() ) );
}

template<typename ValueType>
CUDA_CALLABLE_MEMBER inline Complex<ValueType> fabs( const Complex<ValueType>& a )
{
    Complex<ValueType> x = a;

    if( x.real() < 0 )
    {
        x.real( x.real() * -1 );
    }

    if( x.imag() < 0 )
    {
        x.imag( x.imag() * -1 );
    }

    return x;
}

/**
 * @brief Check if a is lower than b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is lower than b
 */
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, <, int )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, <, long )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, <, float )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<, <, double )
COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator<, <, long double )

/**
 * @brief Check if a is greater than b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is greater than b
 */
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, >, int )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, >, long )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, >, float )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>, >, double )
COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator>, >, long double )

/**
 * @brief Check if a is lower than b or equal to b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is lower than b or equal to b
 */
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <=, int )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <=, long )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <=, float )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator<=, <=, double )
COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator<=, <=, long double )

/**
 * @brief Check if a is greater than b or equal to b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is greater than b or equal to b
 */
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >=, int )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >=, long )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >=, float )
COMPLEX_OPERATOR_COMPARISON_CUDA( operator>=, >=, double )
COMPLEX_OPERATOR_COMPARISON_NONCUDA( operator>=, >=, long double )

/**
 * @brief Check equality of a and b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is equal to b
 */
COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, ==, &&, int )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, ==, &&, long )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, ==, &&, float )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator==, ==, &&, double )
COMPLEX_OPERATOR_EQUALITY_NONCUDA( operator==, ==, &&, long double )

/**
 * @brief Check inequality of a and b.
 *
 * @param[in] a     the 1st Complex to compare this to.
 * @param[in] b     the 2nd Complex to compare this to.
 * @return          if a is unequal to b
 */
COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, !=, ||, int )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, !=, ||, long )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, !=, ||, float )
COMPLEX_OPERATOR_EQUALITY_CUDA( operator!=, !=, ||, double )
COMPLEX_OPERATOR_EQUALITY_NONCUDA( operator!=, !=, ||, long double )

template<typename ValueType>
CUDA_CALLABLE_MEMBER inline Complex<ValueType> operator-( const Complex<ValueType>& a )
{
    Complex<ValueType> x;
    x.real( -a.real() );
    x.imag( -a.imag() );
    return x;
}

/**
 * @brief Add Complex a with Complex b
 *
 * @param[in] a     1st Complex.
 * @param[in] b     2nd Complex.
 * @return          sum
 */
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, +=, int )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, +=, long )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, +=, float )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator+, +=, double )
COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator+, +=, long double )

/**
 * @brief Subtract Complex a with Complex b
 *
 * @param[in] a     1st Complex.
 * @param[in] b     2nd Complex.
 * @return          difference
 */
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -=, int )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -=, long )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -=, float )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator-, -=, double )
COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator-, -=, long double )

/**
 * @brief Multiply Complex a with Complex b
 *
 * @param[in] a     1st Complex.
 * @param[in] b     2nd Complex.
 * @return          product
 */
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *=, int )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *=, long )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *=, float )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator*, *=, double )
COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator*, *=, long double )

/**
 * @brief Divide Complex a with Complex b
 *
 * @param[in] a     1st Complex.
 * @param[in] b     2nd Complex.
 * @return          quotient
 */
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /=, int )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /=, long )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /=, float )
COMPLEX_OPERATOR_NONMEMBER_CUDA( operator/, /=, double )
COMPLEX_OPERATOR_NONMEMBER_NONCUDA( operator/, /=, long double )

template<typename ValueType>
CUDA_CALLABLE_MEMBER inline Complex<ValueType> sqrt( const Complex<ValueType>& a )
{
    ValueType x = a.real();
    ValueType y = a.imag();

    if( x == ValueType() )
    {
        ValueType t = sqrt( abs( y ) / 2 );
        return Complex<ValueType>( t, y < ValueType() ? -t : t );
    }
    else
    {
        ValueType t = sqrt( 2 * ( abs( a ) + abs( x ) ) );
        ValueType u = t / 2;
        return x > ValueType() ? Complex<ValueType>( u, y / t ) :
                        Complex<ValueType>( abs( y ) / t, y < ValueType() ? -u : u );
    }
}

/*
 * long double version must not be CUDA_CALLABLE_MEMBER
 */
inline Complex<long double> sqrt( const Complex<long double>& a )
{
    long double x = a.real();
    long double y = a.imag();

    if( x == 0.0 )
    {
        long double t = sqrt( abs( y ) / 2 );
        return Complex<long double>( t, y < 0.0 ? -t : t );
    }
    else
    {
        long double t = sqrt( 2 * ( abs( a ) + abs( x ) ) );
        long double u = t / 2;
        return x > 0.0 ? Complex<long double>( u, y / t ) : Complex<long double>( abs( y ) / t, y < 0.0 ? -u : u );
    }
}

template<typename ValueType>
std::ostream& operator<<( std::ostream& stream, const Complex<ValueType>& object )
{
    if( object.imag() == 0 )
    {
        stream << object.real();
    }
    else
    {
        stream << object.real() << " " << object.imag();
    }

    return stream;
}

template<typename ValueType,typename InputType1,typename InputType2>
std::basic_istream<InputType1,InputType2>&
operator>>( std::basic_istream<InputType1,InputType2>& input, Complex<ValueType>& x )
{
    ValueType real = 0;
    ValueType imag = 0;
    input >> real;
    input >> imag;
    x.real( real );
    x.imag( imag );
    return input;
}

} /* end namespace lama */

} /* end namespace scai */

#undef CUDA_CALLABLE_MEMBER
#undef COMPLEX_SET_REAL
#undef COMPLEX_SET_COMPLEX
#undef COMPLEX_PLUS_REAL
#undef COMPLEX_PLUS_COMPLEX
#undef COMPLEX_MINUS_REAL
#undef COMPLEX_MINUS_COMPLEX
#undef COMPLEX_MULTIPLY_REAL
#undef COMPLEX_MULTIPLY_COMPLEX
#undef COMPLEX_DIVIDE_REAL
#undef COMPLEX_DIVIDE_COMPLEX
#undef COMPLEX_CAST_REAL
#undef COMPLEX_CAST_COMPLEX
#undef COMPLEX_CONSTRUCTOR_REAL
#undef COMPLEX_CONSTRUCTOR_CUDA
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

