/**
 * @file FileType.hpp
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
 * @brief FileType.hpp
 * @author Kai Buschulte
 * @date 12.05.2010
 * @since 1.0.0
 */

#pragma once

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai::common;

namespace scai
{

namespace lama
{

/** Define an own namespace for enumeration types. */

namespace File
{
/**
 * @brief Defines the supported file types
 */
enum FileType
{
    /**
     * @brief binary format without header informations in the data file
     */
    BINARY,
    /**
     * @brief ascii format
     */
    FORMATTED,
    /**
     * @brief xdr binary format which considers the endianess
     */
    XDR,
    /**
     * @brief the Matrix Market Format
     *        (see http://math.nist.gov/matrixMarket for details).
     */
    MATRIX_MARKET
};

/**
 * @brief Defines the supported data types of the external files
 *
 * For convenience: values and order should be the same as Scalar::ScalarType
 */
enum DataType
{
    INTEGER, //!< synonymous for IndexType
    FLOAT,
    DOUBLE,
    LONG_DOUBLE,
    COMPLEX,
    DOUBLE_COMPLEX,
    LONG_DOUBLE_COMPLEX,
    PATTERN, //!< file for matrix does not contain values
    INTERNAL, //!< takes the internal data type currently used
    UNKNOWN //!< not specified or identified
};

/**
 * @brief Defines the supported index data types of the external files
 */
enum IndexDataType
{
    LONG, INT
};

}
;
// namespace File

/** @brief Help routine to determine the size (in bytes) for the values in a file.
 *
 *  @tparam    ValueType specifies the internal data type take as decision for INTERNAL
 *  @param[in] dataType specifies the file type that is asked for
 *  @returns the size in bytes for the values in the file with the given type
 */

template<typename ValueType>
long getDataTypeSize( const File::DataType dataType )
{
    switch( dataType )
    {
        case File::DOUBLE:
            return TypeTraits<double>::size;

        case File::FLOAT:
            return TypeTraits<float>::size;

        case File::COMPLEX:
            return TypeTraits<ComplexFloat>::size;

        case File::DOUBLE_COMPLEX:
            return TypeTraits<ComplexDouble>::size;

        case File::INTERNAL:
            return TypeTraits<ValueType>::size;

        case File::PATTERN:
            return 0;

        default:
            return -1;
    }
}

/** @brief Determine the file type by its size
 *
 *  Determination of file type by size is ambiguous, e.g. Complex and Double
 *  have same size. If ambiguous, ValueType is the preferred one.
 */

template<typename ValueType>
File::DataType getDataType( const long dataTypeSize )
{
    if( dataTypeSize == TypeTraits<ValueType>::size )
    {
        return File::INTERNAL;
    }
    else if( dataTypeSize == 0 )
    {
        return File::PATTERN;
    }
    else if( dataTypeSize == TypeTraits<float>::size )
    {
        return File::FLOAT;
    }
    else if( dataTypeSize == TypeTraits<double>::size )
    {
        return File::DOUBLE;
    }
    else if( dataTypeSize == TypeTraits<ComplexFloat>::size )
    {
        // never be here as: TypeTraits<double>::size == TypeTraits<ComplexFloat>::size
        // Complex files are only used by INTERNAL
        return File::COMPLEX;
    }
    else if( dataTypeSize == TypeTraits<ComplexDouble>::size )
    {
        return File::DOUBLE_COMPLEX;
    }
    else if( dataTypeSize == TypeTraits<LongDouble>::size )
    {
        // never be here as: TypeTraits<ComplexDouble>::size == TypeTraits<LongDouble>::size
        // LongDouble files are only used by INTERNAL
        return File::LONG_DOUBLE;
    }
    else
    {
        return File::UNKNOWN;
    }
}

static inline
long getIndexDataTypeSize( const File::IndexDataType indexDataType )
{
    switch( indexDataType )
    {
        case File::LONG:
            return TypeTraits<long>::size;

        case File::INT:
            return TypeTraits<int>::size;

        default:
            return 0;
    }
}

} /* end namespace lama */

} /* end namespace scai */
