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
#include <common/config.hpp>

// base classes
#include <common/Printable.hpp>

// others
#include <common/Exception.hpp>

#include <iostream>

namespace common
{

/** Data type that is used for indexing in LAMA.
 *
 *  Currently, it is still necessary that it is a signed data type.
 *  int is the good choice, might be long int or long long int for
 *  future versions that deal with very large matrices even on on processor.
 */
typedef int IndexType;

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

    enum ScalarType
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
    };
}

typedef scalar::ScalarType ScalarType;

/**
 * @brief Conversion of a C type into value of enum ScalarType.
 *
 * @tparam ValueType    C++ type that should be converted
 * @return      value of enum type ScalarType that represents the C++ type.
 */
template<typename ValueType> inline ScalarType getScalarType()
{
    return scalar::UNKNOWN;
}
    
template<>
inline ScalarType getScalarType<IndexType>()
{
    return scalar::INDEX_TYPE;
}

template<>
inline ScalarType getScalarType<float>()
{
    return scalar::FLOAT;
}

template<>
inline ScalarType getScalarType<double>()
{
    return scalar::DOUBLE;
}

} //namespace 

/** Output of ScalarType in stream by writing strings instead of numbers  */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const common::ScalarType& object );
