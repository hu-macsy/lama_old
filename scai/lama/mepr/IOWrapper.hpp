/**
 * @file lama/mepr/IOWrapper.hpp
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
 * @brief Wrapper for IO functions
 * @author Eric Schricker
 * @date 08.03.2016
 */

#pragma once

// local libraries
#include <scai/lama/io/IOUtils.hpp>

// scai internal libraries
#include <scai/common/ScalarType.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace lama {

namespace mepr {

/*
 * Forward declaration
 */
template<typename ValueType, typename TList> struct IOWrapper;

/*
 * Termination
 */
template<typename ValueType>
struct IOWrapper<ValueType, common::mepr::NullType>
{
    static bool readXDR( const common::scalar::ScalarType, XDRFileStream&, ValueType*, const IndexType ) { return false; }
    static bool writeXDR( const common::scalar::ScalarType, XDRFileStream&, const ValueType*, const IndexType ) { return false; }
    static bool readBinary( const common::scalar::ScalarType, std::fstream&, ValueType*, const IndexType ) { return false; }
    static bool writeBinary( const common::scalar::ScalarType, std::fstream&, const ValueType*, const IndexType, const IndexType ) { return false; }

    /*
     *
     */
    static bool readXDR( const long, XDRFileStream&, ValueType*, const IndexType, const IndexType ) { return false; }
    static bool writeXDR( const long, XDRFileStream&, const ValueType*, const IndexType, const IndexType ) { return false; }
    static bool readBinary( const long, std::fstream&, ValueType*, const IndexType, const IndexType ){ return false; }
    static bool writeBinary( const long, std::fstream&, const ValueType*, const IndexType, const IndexType ){ return false; }

};

/*
 * Step n
 */
template<typename ValueType, typename H, typename T>
struct IOWrapper<ValueType, common::mepr::TypeList<H,T> >
{
    static bool readXDR( const common::scalar::ScalarType dataType, XDRFileStream& in, ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::readXDR<H,ValueType>( in, data, n );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::readXDR( dataType, in, data, n );
        }
    }

    static bool writeXDR( const common::scalar::ScalarType dataType, XDRFileStream& out, const ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::writeXDR<H,ValueType>( out, data, n );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::writeXDR( dataType, out, data, n );
        }
    }

    static bool readBinary( const common::scalar::ScalarType dataType, std::fstream& in, ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::readBinaryData<H, ValueType>( in, data, n );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::writeBinary( dataType, in, data, n );
        }
    }

    static bool writeBinary( const common::scalar::ScalarType dataType, std::fstream& out, const ValueType* data, const IndexType n, const IndexType offset = 0 )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::writeBinary<H,ValueType>( out, data, n, offset );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::writeBinary( dataType, out, data, n, offset );
        }
    }

    /*
     *
     */
    static bool readXDR( const long dataTypeSize, XDRFileStream& in, ValueType* data, const IndexType n, const IndexType offset = 0 )
    {
        if( dataTypeSize == sizeof( H ) )
        {
            IOUtils::readXDR<H, ValueType>( in, data, n, offset );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::readXDR( dataTypeSize, in, data, n, offset );
        }
    }

    static bool writeXDR( const long dataTypeSize, XDRFileStream& out, const ValueType* data, const IndexType n, const IndexType offset = 0)
    {
        if( dataTypeSize == sizeof( H ) )
        {
            IOUtils::writeXDR<H, ValueType>( out, data, n, offset );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::writeXDR( dataTypeSize, out, data, n, offset );
        }
    }

    static bool readBinary( const long expectedSize, std::fstream& in, ValueType* data, const IndexType n, const IndexType offset = 0 )
    {
        if( expectedSize == sizeof( H ) )
        {
            IOUtils::readBinaryData<H, ValueType>( in, data, n, offset );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::readBinary( expectedSize, in, data, n, offset );
        }
    }

    static bool writeBinary( const long dataTypeSize, std::fstream& out, const ValueType* data, const IndexType n, const IndexType offset = 0 )
    {
        if( dataTypeSize == sizeof( H ) )
        {
            IOUtils::writeBinary<H, ValueType>( out, data, n, offset );
            return true;
        }
        else
        {
            return IOWrapper<ValueType, T>::writeBinary( dataTypeSize, out, data, n, offset );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
