/**
 * @file ScalarType.hpp
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
 * @brief Definition of enum type for scalar value types typically used in numerical applications.
 * @author Thomas Brandes
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/common/macros/throw.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
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
enum class ScalarType
{
    INT,                 //!<  synonymous for int (signed 32 bit)
    LONG,                //!<  synonymous for long (signed 64 bit)
    FLOAT,               //!<  synonymous for float
    DOUBLE,              //!<  synonymous for double
    LONG_DOUBLE,         //!<  synonymous for long double
    COMPLEX,             //!<  synonymous for complex
    DOUBLE_COMPLEX,      //!<  synonymous for double complex
    LONG_DOUBLE_COMPLEX, //!<  synonymous for long double complex
    CHAR,                //!<  synonymous for char(signed 8 bit)
    UNSIGNED_INT,        //!<  synonymous for uint(unsigned 32 bit)
    UNSIGNED_LONG,       //!<  synonymous for size_t(unsigned 64 bit)

    INDEX_TYPE,          //!<  take the type as defined for IndexType, getScalarType<IndexType>()
    PATTERN,             //!<  dummy type of size 0
    INTERNAL,            //!<  take the type currently in use, getScalarType<ValueType>()
    UNKNOWN
}; 

/** Simple enun class for the two parts of a complex value */

enum class ComplexPart
{
    REAL,   // real part of a complex number
    IMAG    // imaginary part of a complex number
};

/** Query an enum value whether it is complex or not. */

COMMON_DLL_IMPORTEXPORT bool isComplex( const ScalarType stype );

/** Query an enum value whether it is supported and numeric. */

COMMON_DLL_IMPORTEXPORT bool isNumeric( const ScalarType stype );

/** Query an enum value about the precision used in output */

COMMON_DLL_IMPORTEXPORT int precision( const ScalarType stype );

/** Query sizeof( ScalarType stype ) same as sizeof( type ) for which stype stands for */

COMMON_DLL_IMPORTEXPORT size_t typeSize( const ScalarType stype );

/*
 * Output of ScalarType in stream by writing strings instead of numbers
 */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ScalarType& object );

/**
 * @brief Conversion of a ScalarType to string
 *
 */

COMMON_DLL_IMPORTEXPORT const char* scalarType2str( const ScalarType stype );

COMMON_DLL_IMPORTEXPORT ScalarType str2ScalarType( const char* str );

} /* end namespace common */

} /* end namespace scai */
