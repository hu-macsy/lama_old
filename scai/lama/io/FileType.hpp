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
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>

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
    MATRIX_MARKET, 
    /**
     * @brief unspecified, used internally 
     */
    DEFAULT
};

/*
 * Output of ScalarType in stream by writing strings instead of numbers
 */

static inline std::ostream& operator<<( std::ostream& stream, const FileType& object )
{
    switch ( object ) 
    {
        case BINARY:
            stream << "BINARY";
            break;
        case FORMATTED:
            stream << "FORMATTED";
            break;
        case XDR:
            stream << "XDR";
            break;
        case MATRIX_MARKET:
            stream << "MATRIX_MARKET";
            break;
        default:
            stream << "<unknown_file_type>";
     }
     return stream;
}

// typedef ScalarType DataType;

/**
 * @brief Defines the supported index data types of the external files
 */
enum IndexDataType
{
    LONG, INT
};

}  // namespace File

/** @brief Help routine to determine the size (in bytes) for the values in a file.
 *
 *  @tparam    ValueType specifies the internal data type take as decision for INTERNAL
 *  @param[in] dataType specifies the file type that is asked for
 *  @returns the size in bytes for the values in the file with the given type
 */

template<typename ValueType>
long getDataTypeSize( const common::scalar::ScalarType dataType )
{
    long s = common::mepr::ScalarTypeHelper<ARITHMETIC_HOST_LIST>::sizeOf( dataType );

    if( s != 0 )
    {
        return s;
    }

    switch( dataType )
    { 
        case common::scalar::INTERNAL:
            return sizeof( ValueType );

        case common::scalar::PATTERN:
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
common::scalar::ScalarType getDataType( const long dataTypeSize )
{
    if( dataTypeSize == sizeof( ValueType ) )
    {
        return common::scalar::INTERNAL;
    }
    else if( dataTypeSize == 0 )
    {
        return common::scalar::PATTERN;
    }

    return common::mepr::ScalarTypeHelper<ARITHMETIC_HOST_LIST>::getBySize( dataTypeSize );
}

static inline
long getIndexDataTypeSize( const File::IndexDataType indexDataType )
{
    switch( indexDataType )
    {
        case File::LONG:
            return sizeof( long );

        case File::INT:
            return sizeof( int );

        default:
            return 0;
    }
}

} /* end namespace lama */

} /* end namespace scai */
