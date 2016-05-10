/**
 * @file XDRFileStream.hpp
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
 * @brief XDRFileStream.hpp
 * @author Kai Buschulte
 * @date 15.04.2010
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>

// std
#include <fstream>
#include <cstdlib>
#include <typeinfo>
#include <sstream>

namespace scai
{

namespace lama
{

/**
 * @brief XDR file stream class to access XDR files
 * Supported data types: (unsigned) int, (unsigned) long, float, double
 */
class COMMON_DLL_IMPORTEXPORT XDRFileStream
{
public:
    /**
     * @brief   The default constructor creates a XDR stream in the specified
     *          openmode. (in/out)
     *
     * @param[in] filename  Includes (the path and) the file name.
     * @param[in] openmode  Openmode of the XDR stream defined in standard C++
     *                  I/O openmode (std::ios::in or std::ios::out). The binary
     *                  (std::ios::in) mode is automatically added.
     */
    XDRFileStream( const std::string& filename, const std::ios_base::openmode openmode );

    /**
     * @brief The default destructor destroys the XDR stream if the file has
     *        been opened.
     */
    virtual ~XDRFileStream();

    /**
     * @brief Reads Values(s) of Type ValueType from the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be read from
     *                     the XDR stream.
     * @param[in]  n       Number of values which should be read to the input field.
     */
    template<typename ValueType>
    void read( ValueType * const input, const std::streamsize n );

    /**
     * @brief Reads Values(s) of Type ValueType from the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be read from
     *                     the XDR stream.
     */
    template<typename ValueType>
    void read( ValueType * const input );

    /**
     * @brief Writes Values of Type ValueType to the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be written to the
     *                     XDR stream.
     * @param[in] n        Number of values which should be written to the output
     */
    template<typename ValueType>
    void write( ValueType * const input, const std::streamsize n );

    /**
     * @brief Writes Values of Type ValueType to the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be written to the
     *                     XDR stream.
     */
    template<typename ValueType>
    void write( ValueType * const input );

    /**
     * @brief Checks if the vector file could be opened.
     *
     * @return Returns true if the file could be successfully opened.
     */
    bool is_open();

    /**
     * @brief Destroys the XDR stream if it has been opened.
     */
    void close();

private:

    std::string m_filename;
    std::fstream m_filestream;

    std::ios_base::openmode m_openmode;

    void openFilePointer( const std::string filename, const std::ios_base::openmode openmode );
    template<typename ValueType>
    void xdrRead( ValueType* const data );
    template<typename ValueType>
    void xdrWrite( const ValueType* const data );

    template<typename ValueType>
    int getSize( const ValueType );

    bool isLittleEndian();

};
//class XDRFileStream

template<typename ValueType>
void XDRFileStream::write( ValueType * const input, const std::streamsize n )
{
    if( m_openmode & std::ios::in )
    {
        COMMON_THROWEXCEPTION( "XDRFileStream: Stream is not in Output mode" )
    }

    for( std::streamsize i = 0; i < n; i++ )
    {
        xdrWrite( &input[i] );
    }
}

template<typename ValueType>
void XDRFileStream::write( ValueType * const input )
{
    write( input, 1 );
}

template<typename ValueType>
void XDRFileStream::read( ValueType * const input, const std::streamsize n )
{
    if( m_openmode & std::ios::out )
    {
        COMMON_THROWEXCEPTION( "XDRFileStream: Stream is not in Input mode" )
    }

    for( std::streamsize i = 0; i < n; i++ )
    {
        xdrRead( &input[i] );
    }
}

template<typename ValueType>
void XDRFileStream::read( ValueType * const input )
{
    read( input, 1 );
}

template<typename ValueType>
void XDRFileStream::xdrRead( ValueType* const data )
{
    if( !is_open() )
    {
        COMMON_THROWEXCEPTION( "Error trying to read from closed XDRFileStream" )
    }

    const int length = getSize( *data );
    *data = 0;
    char* dataptr = reinterpret_cast<char*>( data );

    if( isLittleEndian() )
    {
        for( int pos = length - 1; pos >= 0; pos-- )
        {
            m_filestream.read( dataptr + pos, 1 );
        }
    }
    else
    {
        m_filestream.read( dataptr, length );
    }

    //correct the signed bit
    if( typeid(ValueType) == typeid(long) )
    {
        int i = 0 | static_cast<long>( *data );
        *data = 0;
        *data = static_cast<ValueType>( i );
    }
}
template<typename ValueType>
void XDRFileStream::xdrWrite( const ValueType* const data )
{
    if( !is_open() )
    {
        COMMON_THROWEXCEPTION( "Error trying to read from closed XDRFileStream" )
    }

    const int length = getSize( *data );
    //correct the signed bit
    ValueType tempData = *data;

    if( typeid(ValueType) == typeid(unsigned long) )
    {
        //Check if unsigned long is over the limit of XDR
        if( ( ( static_cast<long>( *data ) ) ) > ( ( static_cast<long>( 1 ) << ( length * 8 ) ) - 1 ) )
        {
            COMMON_THROWEXCEPTION( "unsigned long is to big for XDR (Limit 4 Byte)" )
        }
    }

    if( typeid(ValueType) == typeid(long) )
    {
        if( std::abs( static_cast<long>( *data ) ) >= ( 1 << 30 ) - 1 )
        {
            COMMON_THROWEXCEPTION( "long is to big for XDR (Limit 4 Byte)" )
        }

        int temp = static_cast<long>( *data );
        tempData = static_cast<ValueType>( 0 | temp );
    }

    char* dataptr = reinterpret_cast<char*>( &tempData );

    if( isLittleEndian() )
    {
        for( int pos = length - 1; pos >= 0; pos-- )
        {
            m_filestream.write( dataptr + pos, 1 );
        }
    }
    else
    {
        m_filestream.write( dataptr, length );
    }
}

template<typename ValueType>
int XDRFileStream::getSize( const ValueType )
{
    int size = 0;

    if( typeid(ValueType) == typeid(double) )
    {
        size = 8;
    }
    else if( typeid(ValueType) == typeid(float) )
    {
        size = 4;
    }
    else if( typeid(ValueType) == typeid(int) )
    {
        size = 4;
    }
    else if( typeid(ValueType) == typeid(long) )
    {
        size = 4;
    }
    else if( typeid(ValueType) == typeid(unsigned long) )
    {
        size = 4;
    }
    else if( typeid(ValueType) == typeid(unsigned int) )
    {
        size = 4;
    }
    else
    {
        COMMON_THROWEXCEPTION( "XDRFileStream: Type not permitted: " << typeid(ValueType).name() )
    }

    return size;
}

} /* end namespace lama */

} /* end namespace scai */
