/**
 * @file lama/mepr/IOWrapper.hpp
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
 * @brief Wrapper for IO functions
 * @author Eric Schricker
 * @date 08.03.2016
 */

#pragma once

#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/io/FileIO.hpp>

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
            FileIO::readBinaryData<H, ValueType>( in, data, n );
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
            FileIO::readBinaryData<H, ValueType>( in, data, n, offset );
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
