/**
 * @file FileType.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief FileType.hpp
 * @author Kai Buschulte
 * @date 12.05.2010
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
    long s = common::mepr::ScalarTypeHelper<SCAI_ARITHMETIC_HOST_LIST>::sizeOf( dataType );

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

    return common::mepr::ScalarTypeHelper<SCAI_ARITHMETIC_HOST_LIST>::getBySize( dataTypeSize );
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
