/*
 * @file Complex.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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

#ifndef COMPLEX_H_
#define COMPLEX_H_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__ __host__
#include <cuComplex.h>
#include <math.h>
#else
#define CUDA_CALLABLE_MEMBER
#include <cmath>
#endif

#include <lama/exception/LAMAAssert.hpp>

/*
 * Macros used for building the operators, providing the functionality
 * Every macro is defined with the suffix REAL and COMPLEX. Each is
 * providing another functionality.
 */

#define COMPLEX_SET_REAL                                                                            \
real( static_cast<U>( t ) );                                                                        \
imag( static_cast<U>( 0 ) );                                                                        \
return *this;

#define COMPLEX_SET_COMPLEX                                                                         \
real( static_cast<U>( t.real() ) );                                                                 \
imag( static_cast<U>( t.imag() ) );                                                                 \
return *this;

#define COMPLEX_PLUS_REAL                                                                           \
real( real() + static_cast<U>( t ) );                                                               \
return *this;

#define COMPLEX_PLUS_COMPLEX                                                                        \
real( real() + static_cast<U>( t.real() ) );                                                        \
imag( imag() + static_cast<U>( t.imag() ) );                                                        \
return *this;

#define COMPLEX_MINUS_REAL                                                                          \
real( real() - static_cast<U>( t ) );                                                               \
return *this;

#define COMPLEX_MINUS_COMPLEX                                                                       \
real( real() - static_cast<U>( t.real() ) );                                                        \
imag( imag() - static_cast<U>( t.imag() ) );                                                        \
return *this;

#define COMPLEX_MULTIPLY_REAL                                                                       \
real( real() * static_cast<U>( t ) );                                                               \
imag( imag() * static_cast<U>( t ) );                                                               \
return *this;

#define COMPLEX_MULTIPLY_COMPLEX                                                                    \
U a = real();                                                                                       \
U b = imag();                                                                                       \
U c = static_cast<U>( t.real() );                                                                   \
U d = static_cast<U>( t.imag() );                                                                   \
real( a * c - b * d );                                                                              \
imag( a * d + b * c );                                                                              \
return *this;

#define COMPLEX_DIVIDE_REAL                                                                         \
real( real() / static_cast<U>( t ) );                                                               \
imag( imag() / static_cast<U>( t ) );                                                               \
return *this;

#define COMPLEX_DIVIDE_COMPLEX                                                                      \
U a = real();                                                                                       \
U b = imag();                                                                                       \
U c = static_cast<U>( t.real() );                                                                   \
U d = static_cast<U>( t.imag() );                                                                   \
U c2d2 = ( c * c + d * d);                                                                          \
real( ( a * c + b * d ) / ( c2d2 ) );                                                               \
imag( ( b * c - a * d ) / ( c2d2 ) );                                                               \
return *this;

#define COMPLEX_CAST_REAL(type)                                                                     \
return static_cast<type>( metrikCuda() );

#define COMPLEX_CAST_COMPLEX(type)                                                                  \
return Complex<type>( static_cast<type>( real() ), static_cast<type>( imag() ) );

#define COMPLEX_CONSTRUCTOR_REAL                                                                    \
this->real( static_cast<U>( value ) );                                                              \
this->imag( static_cast<U>( 0 ) );

/*
 * Macros which build the operators, the macros above will be used for the functionality
 * The macros with the suffix CUDA are getting in an cuda environment the prefix __host__ __device__
 * The macros with the suffix NONCUDA are getting no prefix
 *
 * For creating constructors
 *
 */

#define COMPLEX_CONSTRUCTOR_CUDA( type, method )                                                    \
CUDA_CALLABLE_MEMBER                                                                                \
inline Complex( const type value )                                                                  \
{                                                                                                   \
    method                                                                                          \
}

#define COMPLEX_CONSTRUCTOR_NONCUDA( type, method )                                                 \
		                                                                                            \
inline Complex( const type value )                                                                  \
{                                                                                                   \
    method                                                                                          \
}

/*
 * For member calculation operators =, +=, -=, *= and /=
 */

#define COMPLEX_OPERATOR_CUDA( op, type, method )                                                   \
CUDA_CALLABLE_MEMBER inline Complex<U>& op( const type t )                                          \
{                                                                                                   \
    method                                                                                          \
}                                                                                                   \
CUDA_CALLABLE_MEMBER inline volatile Complex<U>& op( const type t ) volatile                        \
{                                                                                                   \
    method                                                                                          \
}

#define COMPLEX_OPERATOR_NONCUDA( op, type, method )                                                \
inline Complex<U>& op( const type t )                                                               \
{                                                                                                   \
    method                                                                                          \
}                                                                                                   \
inline volatile Complex<U>& op( const type t ) volatile                                             \
{                                                                                                   \
    method                                                                                          \
}

/*
 * For cast operators
 */

#define COMPLEX_OPERATOR_CAST_CUDA( type, method )                                                  \
CUDA_CALLABLE_MEMBER operator type() const                                                          \
{                                                                                                   \
    method                                                                                          \
}

#define COMPLEX_OPERATOR_CAST_NONCUDA( type, method )                                               \
operator type() const                                                                               \
{                                                                                                   \
    method                                                                                          \
}

/*
 * For comparison operators: <, >, <=, >=
 *
 */

#define COMPLEX_OPERATOR_COMPARISON_CUDA(op, sign, type)                                            \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const Complex<type>& b )                  \
{                                                                                                   \
    return a.metrikCuda() sign static_cast<T>(b.metrikCuda());                                      \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const type& b )                           \
{                                                                                                   \
    return a.metrikCuda() sign static_cast<T>(b);                                                   \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<T>& b )                           \
{                                                                                                   \
    return a sign static_cast<type>(b.metrikCuda());                                                \
}

#define COMPLEX_OPERATOR_COMPARISON_NONCUDA(op, sign, type)                                         \
template<typename T>                                                                                \
inline bool op( const Complex<T>& a, const Complex<type>& b )                                       \
{                                                                                                   \
    return a.metrikHost() sign static_cast<T>(b.metrikHost());                                      \
}                                                                                                   \
template<typename T>                                                                                \
inline bool op( const Complex<T>& a, const type& b )                                                \
{                                                                                                   \
    return a.metrikHost() sign static_cast<T>(b);                                                   \
}                                                                                                   \
template<typename T>                                                                                \
inline bool op( const type& a, const Complex<T>& b )                                                \
{                                                                                                   \
    return a sign static_cast<type>(b.metrikHost());                                                \
}

/*
 * For equality operators: ==, !=
 */

#define COMPLEX_OPERATOR_EQUALITY_CUDA(op, sign, connection, type)                                  \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const Complex<type>& b )                  \
{                                                                                                   \
    return a.real() sign b.real() connection a.imag() sign b.imag();                                \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const type& b )                           \
{                                                                                                   \
    return a.real() sign b connection a.imag() sign 0;                                              \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<T>& b )                           \
{                                                                                                   \
    return a sign b.real() connection 0 sign b.imag();                                              \
}

#define COMPLEX_OPERATOR_EQUALITY_NONCUDA(op, sign, connection, type)                               \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const Complex<type>& b )                  \
{                                                                                                   \
    return a.real() sign b.real() connection a.imag() sign b.imag();                                \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const Complex<T>& a, const type& b )                           \
{                                                                                                   \
    return a.real() sign b connection a.imag() sign 0;                                              \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline bool op( const type& a, const Complex<T>& b )                           \
{                                                                                                   \
    return a sign b.real() connection 0 sign b.imag();                                              \
}

/*
 * For non-member calculation operators: +,-,*,/
 */

#define COMPLEX_OPERATOR_NONMEMBER_CUDA(op, sign, type)                                             \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline Complex<T> op( const Complex<T>& a, const Complex<type>& b )            \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline Complex<T> op( volatile Complex<T>& a, const Complex<type>& b )         \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
CUDA_CALLABLE_MEMBER inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )   \
{                                                                                                   \
    Complex<type> x = a;                                                                            \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline Complex<T> op( const Complex<T>& a, const type& b )                     \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign static_cast<T>(b);                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline Complex<T> op( const type& a, const Complex<T>& b )                     \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
CUDA_CALLABLE_MEMBER inline Complex<T> op( Complex<T>& a, Complex<type>& b )                        \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \

#define COMPLEX_OPERATOR_NONMEMBER_NONCUDA(op, sign, type)                                          \
template<typename T>                                                                                \
inline Complex<T> op( const Complex<T>& a, const Complex<type>& b )                                 \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
inline Complex<T> op( volatile Complex<T>& a, const Complex<type>& b )                              \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
inline Complex<type> op( volatile Complex<type>& a, const Complex<type>& b )                        \
{                                                                                                   \
    Complex<type> x = a;                                                                            \
    x sign b;                                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
inline Complex<T> op( const Complex<T>& a, const type& b )                                          \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
    x sign static_cast<T>(b);                                                                       \
    return x;                                                                                       \
}                                                                                                   \
template<typename T>                                                                                \
inline Complex<T> op( const type& a, const Complex<T>& b )                                          \
{                                                                                                   \
    Complex<T> x = a;                                                                               \
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

template<typename U>
class LAMA_DLL_IMPORTEXPORT Complex
{
public:
    CUDA_CALLABLE_MEMBER
    explicit Complex();

    template<typename T>
    CUDA_CALLABLE_MEMBER explicit inline Complex( const T value );

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
    template<typename S,typename T>
    CUDA_CALLABLE_MEMBER inline Complex( const S real, const T imag )
    {
        this->real( static_cast<U>( real ) );
        this->imag( static_cast<U>( imag ) );
    }

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @param[in]   value the value this complex should represent
     */
    template<typename T>
    CUDA_CALLABLE_MEMBER inline Complex( const Complex<T>& value )
    {
        this->real( static_cast<U>( value.real() ) );
        this->imag( static_cast<U>( value.imag() ) );
    }

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @param[in]   value the value this complex should represent
     */
    CUDA_CALLABLE_MEMBER
    inline Complex( volatile const Complex<U>& value )
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

    template<typename T>
    CUDA_CALLABLE_MEMBER
    operator Complex<T>() const
    {
        return Complex<T>(static_cast<T>(real()), static_cast<T>(imag()));
    }

    /**
     * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is callable from CUDA space.
     *
     * @return      the metric of a complex number
     */
    CUDA_CALLABLE_MEMBER U metrikCuda( void ) const;

    /**
     * @brief Returns the metric for a complex number which is used to determine an order for complex numbers. Through this it is possible to use functions like min or max. This one is not callable from CUDA space and used for long double operations.
     *
     * @return      the metric of a complex number
     */
    U metrikHost( void ) const;

    /**
     * @brief Returns imaginary part of the complex number
     *
     * @return     the imaginary part of a complex number
     */
    CUDA_CALLABLE_MEMBER
    inline U imag( void )
    {
        return mParts[1];
    }

    CUDA_CALLABLE_MEMBER
    inline U imag( void ) const
    {
        return mParts[1];
    }

    CUDA_CALLABLE_MEMBER
    inline U imag( void ) volatile
    {
        return mParts[1];
    }

    CUDA_CALLABLE_MEMBER
    inline U imag( void ) const volatile
    {
        return mParts[1];
    }

    /**
     * @brief Sets imaginary part of the complex number
     *
     * @param[in]     the new imaginary part of the complex number
     */
    CUDA_CALLABLE_MEMBER
    inline void imag( const U a )
    {
        mParts[1] = a;
    }

    CUDA_CALLABLE_MEMBER
    inline void imag( const U a ) volatile
    {
        mParts[1] = a;
    }

    /**
     * @brief Returns real part of the complex number
     *
     * @return     the real part of a complex number
     */
    CUDA_CALLABLE_MEMBER
    inline U real( void )
    {
        return mParts[0];
    }

    CUDA_CALLABLE_MEMBER
    inline U real( void ) const
    {
        return mParts[0];
    }

    CUDA_CALLABLE_MEMBER
    inline U real( void ) volatile
    {
        return mParts[0];
    }

    CUDA_CALLABLE_MEMBER
    inline U real( void ) const volatile
    {
        return mParts[0];
    }

    /**
     * @brief Sets real part of the complex number
     *
     * @param[in]     the new real part of the complex number
     */
    CUDA_CALLABLE_MEMBER
    inline void real( const U a )
    {
        mParts[0] = a;
    }

    CUDA_CALLABLE_MEMBER
    inline void real( const U a ) volatile
    {
        mParts[0] = a;
    }

private:
    U mParts[2]; // mParts[0] <-- Real, //mParts[1] <-- Imag
};

template<typename T>
CUDA_CALLABLE_MEMBER Complex<T>::Complex()
{
    //Constructor have to be empty. Otherwise cuda will output a lot of warnings
}

template<typename T>
CUDA_CALLABLE_MEMBER T Complex<T>::metrikCuda( void ) const
{
    return sqrt( real() * real() + imag() * imag() );
}

template<typename T>
T Complex<T>::metrikHost( void ) const
{
    return sqrt( real() * real() + imag() * imag() );
}

template<typename T>
CUDA_CALLABLE_MEMBER inline T abs( const Complex<T>& a )
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

template<typename T>
CUDA_CALLABLE_MEMBER inline Complex<T> fabs( const Complex<T>& a )
{
    Complex<T> x = a;
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

template<typename T>
CUDA_CALLABLE_MEMBER inline Complex<T> operator-( const Complex<T>& a )
{
    Complex<T> x;
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

template<typename T>
CUDA_CALLABLE_MEMBER inline Complex<T> sqrt( const Complex<T>& a )
{
    T x = a.real();
    T y = a.imag();
    if( x == T() )
    {
        T t = sqrt( abs( y ) / 2 );
        return Complex<T>( t, y < T() ? -t : t );
    }
    else
    {
        T t = sqrt( 2 * ( abs( a ) + abs( x ) ) );
        T u = t / 2;
        return x > T() ? Complex<T>( u, y / t ) : Complex<T>( abs( y ) / t, y < T() ? -u : u );
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

template<typename T>
std::ostream& operator<<( std::ostream& stream, const Complex<T>& object )
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

template<typename T,typename U,typename V>
std::basic_istream<U,V>&
operator>>( std::basic_istream<U,V>& input, Complex<T>& x )
{
    T real = 0;
    T imag = 0;
    input >> real;
    input >> imag;
    x.real( real );
    x.imag( imag );
    return input;

}

} //namespace lama

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

#endif /* COMPLEX_H_ */