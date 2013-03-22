/**
 * @file XDRFileStream.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief XDRFileStream.hpp
 * @author kbuschul
 * @date 15.04.2010
 * $Id$
 */

#ifndef XDRFILESTREAM_HPP_
#define XDRFILESTREAM_HPP_

// for dll_import
#include <lama/config.hpp>
#include <lama/exception/Exception.hpp>

// others
#include <fstream>
#include <cstdlib>
#include <typeinfo>
#include <sstream>

namespace lama
{

/**
 * @brief XDR file stream class to access XDR files
 * Supported data types: (unsigned) int, (unsigned) long, float, double
 */
class LAMA_DLL_IMPORTEXPORT XDRFileStream
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
     * @brief Reads Values(s) of Type T from the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be read from
     *                     the XDR stream.
     * @param[in]  n       Number of values which should be read to the input field.
     */
    template<typename T>
    void read( T * const input, const std::streamsize n );

    /**
     * @brief Reads Values(s) of Type T from the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be read from
     *                     the XDR stream.
     */
    template<typename T>
    void read( T * const input );

    /**
     * @brief Writes Values of Type T to the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be written to the
     *                     XDR stream.
     * @param[in] n        Number of values which should be written to the output
     */
    template<typename T>
    void write( T * const input, const std::streamsize n );

    /**
     * @brief Writes Values of Type T to the XDR stream.
     *
     * @param[in] input    Data pointer for the data which should be written to the
     *                     XDR stream.
     */
    template<typename T>
    void write( T * const input );

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
    template<typename T>
    void xdrRead( T* const data );
    template<typename T>
    void xdrWrite( const T* const data );

    template<typename T>
    int getSize( const T );

    bool isLittleEndian();

};
//class XDRFileStream

template<typename T>
void XDRFileStream::write( T * const input, const std::streamsize n )
{
    if ( m_openmode & std::ios::in )
    {
        LAMA_THROWEXCEPTION( "XDRFileStream: Stream is not in Output mode" )
    }
    for ( std::streamsize i = 0; i < n; i++ )
    {
        xdrWrite( &input[i] );
    }
}

template<typename T>
void XDRFileStream::write( T * const input )
{
    write( input, 1 );
}

template<typename T>
void XDRFileStream::read( T * const input, const std::streamsize n )
{
    if ( m_openmode & std::ios::out )
    {
        LAMA_THROWEXCEPTION( "XDRFileStream: Stream is not in Input mode" )
    }
    for ( std::streamsize i = 0; i < n; i++ )
    {
        xdrRead( &input[i] );
    }
}

template<typename T>
void XDRFileStream::read( T * const input )
{
    read( input, 1 );
}

template<typename T>
void XDRFileStream::xdrRead( T* const data )
{
    if ( !is_open() )
    {
        LAMA_THROWEXCEPTION( "Error trying to read from closed XDRFileStream" )
    }
    const int length = getSize( *data );
    *data = 0;
    char* dataptr = reinterpret_cast<char*>( data );
    if ( isLittleEndian() )
    {
        for ( int pos = length - 1; pos >= 0; pos-- )
        {
            m_filestream.read( dataptr + pos, 1 );
        }
    }
    else
    {
        m_filestream.read( dataptr, length );
    }
    //correct the signed bit
    if ( typeid(T) == typeid(long) )
    {
        int i = 0 | static_cast<long>( *data );
        *data = 0;
        *data = static_cast<T>( i );
    }
}
template<typename T>
void XDRFileStream::xdrWrite( const T* const data )
{
    if ( !is_open() )
    {
        LAMA_THROWEXCEPTION( "Error trying to read from closed XDRFileStream" )
    }
    const int length = getSize( *data );
    //correct the signed bit
    T tempData = *data;
    if ( typeid(T) == typeid(unsigned long) )
    {
        //Check if unsigned long is over the limit of XDR
        if ( ( ( static_cast<long>( *data ) ) ) > ( ( static_cast<long>( 1 ) << ( length * 8 ) ) - 1 ) )
        {
            LAMA_THROWEXCEPTION( "unsigned long is to big for XDR (Limit 4 Byte)" )
        }
    }
    if ( typeid(T) == typeid(long) )
    {
        if ( std::abs( static_cast<long>( *data ) ) >= ( 1 << 30 ) - 1 )
        {
            LAMA_THROWEXCEPTION( "long is to big for XDR (Limit 4 Byte)" )
        }
        int temp = static_cast<long>( *data );
        tempData = static_cast<T>( 0 | temp );
    }
    char* dataptr = reinterpret_cast<char*>( &tempData );
    if ( isLittleEndian() )
    {
        for ( int pos = length - 1; pos >= 0; pos-- )
        {
            m_filestream.write( dataptr + pos, 1 );
        }
    }
    else
    {
        m_filestream.write( dataptr, length );
    }
}

template<typename T>
int XDRFileStream::getSize( const T )
{
    int size = 0;

    if ( typeid(T) == typeid(double) )
    {
        size = 8;
    }
    else if ( typeid(T) == typeid(float) )
    {
        size = 4;
    }
    else if ( typeid(T) == typeid(int) )
    {
        size = 4;
    }
    else if ( typeid(T) == typeid(long) )
    {
        size = 4;
    }
    else if ( typeid(T) == typeid(unsigned long) )
    {
        size = 4;
    }
    else if ( typeid(T) == typeid(unsigned int) )
    {
        size = 4;
    }
    else
    {
        LAMA_THROWEXCEPTION( (std::string("XDRFileStream: Type not permitted: ") + typeid(T).name()).c_str() )
    }
    return size;
}

} //namespace lama

#endif // XDRFILESTREAM_HPP_
