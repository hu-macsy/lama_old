/**
 * @file TypeTraits.hpp
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
 * @brief Contains the template TypeTraits and specializations.
 * @author Jiri Kraus
 * @date 14.04.2010
 * @since 1.0.0
 */

#ifndef TYPETRAITS_HPP_
#define TYPETRAITS_HPP_

namespace lama
{

/**
 * @brief The template class TypeTraits determines the index and value types.
 *
 * The template class TypeTraits is used to determine the type of indices and
 * values in higher level types like CSRSparseMatrix. Its facilities are used by
 * the template function maxNorm (@see MaxNorm.hpp) among others.
 *
 * @tparam T The type for values and indices.
 */
template<typename T>
class TypeTraits
{
public:
    /**
     * @brief The type for indices.
     */
    typedef typename T::IndexType IndexType;
    /**
     * @brief The type for values.
     */
    typedef typename T::ValueType ValueType;

    typedef typename T::ExpressionMemberType ExpressionMemberType;

    static const long size = sizeof(T);
};

template<>
class TypeTraits<long>
{
public:
    typedef long IndexType;
    typedef long ValueType;

    typedef const long ExpressionMemberType;

    static const long size = 8;
};

template<>
class TypeTraits<int>
{
public:
    typedef int IndexType;
    typedef int ValueType;
    typedef const int ExpressionMemberType;
    static const long size = 4;
};

template<>
class TypeTraits<float>
{
public:
    typedef int IndexType;
    typedef float ValueType;
    typedef const float ExpressionMemberType;
    static const long size = sizeof(float);
};

template<>
class TypeTraits<double>
{
public:
    typedef int IndexType;
    typedef double ValueType;
    typedef const double ExpressionMemberType;
    static const long size = sizeof(double);
};

template<>
class TypeTraits<ComplexFloat>
{
public:
    typedef int IndexType;
    typedef ComplexFloat ValueType;
    typedef const ComplexFloat ExpressionMemberType;
    static const long size = sizeof(ComplexFloat);
};

template<>
class TypeTraits<ComplexDouble>
{
public:
    typedef int IndexType;
    typedef ComplexDouble ValueType;
    typedef const ComplexDouble ExpressionMemberType;
    static const long size = sizeof(ComplexDouble);
};

template<>
class TypeTraits<LongDouble>
{
public:
    typedef int IndexType;
    typedef LongDouble ValueType;
    typedef const LongDouble ExpressionMemberType;
    static const long size = sizeof(LongDouble);
};

} //namespace lama

#endif // TYPETRAITS_HPP_
