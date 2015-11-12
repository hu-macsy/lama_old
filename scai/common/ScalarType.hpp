/**
 * @file ScalarType.hpp
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
 * @brief Definition of enum type for scalar value types typically used in numerical applications.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/common/SCAITypes.hpp>

#include <scai/common/macros/throw.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
{

/** Namespace for enumeration type. */

namespace scalar
{
    /** Enumeration type for supported value types in LAMA.
     *
     *  This enumeration type is result of many query operations for LAMA classes
     *  and avoids expensive calls of the typeid operator.
     *
     *  \code
     *    CSRSparseMatrix<double> a;
     *    a.getValueType()  // returns ScalarType::DOUBLE
     *    getScalarType<double>()  // return ScalarType DOUBLE
     *  \endcode
     *
     *  It is especially useful when casting variables of base classes to derived classes.
     */

    typedef enum
    {
        INDEX_TYPE, //!<  synonymous for IndexType
        FLOAT, //!<  synonymous for float
        DOUBLE, //!<  synonymous for double
        LONG_DOUBLE, //!<  synonymous for long double
        COMPLEX, //!<  synonymous for complex
        DOUBLE_COMPLEX, //!<  synonymous for double complex
        LONG_DOUBLE_COMPLEX, //!<  synonymous for long double complex
        PATTERN, //!<  dummy type of size 0
        INTERNAL, //!<  take the type currently in use, getScalarType<ValueType>()
        UNKNOWN
    } ScalarType;

    /*
     * Output of ScalarType in stream by writing strings instead of numbers
     */
    
    COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ScalarType& object );

} /* end namespace scalar */

// make ScalarType visible, but not the enum values

using scalar::ScalarType;

/**
 * @brief Conversion of a C type into value of enum ScalarType.
 *
 * @tparam ValueType    C++ type that should be converted
 * @return      value of enum type ScalarType that represents the C++ type.
 */
template<typename ValueType> inline scalar::ScalarType getScalarType()
{
    return scalar::UNKNOWN;
}
    
template<>
inline scalar::ScalarType getScalarType<IndexType>()
{
    return scalar::INDEX_TYPE;
}

template<>
inline scalar::ScalarType getScalarType<float>()
{
    return scalar::FLOAT;
}

template<>
inline scalar::ScalarType getScalarType<double>()
{
    return scalar::DOUBLE;
}

template<>
inline scalar::ScalarType getScalarType<LongDouble>()
{
    return scalar::LONG_DOUBLE;
}

template<>
inline scalar::ScalarType getScalarType<ComplexFloat>()
{
    return scalar::COMPLEX;
}

template<>
inline scalar::ScalarType getScalarType<ComplexDouble>()
{
    return scalar::DOUBLE_COMPLEX;
}

template<>
inline scalar::ScalarType getScalarType<ComplexLongDouble>()
{
    return scalar::LONG_DOUBLE_COMPLEX;
}

} /* end namespace common */

} /* end namespace scai */
