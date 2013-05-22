/**
 * @file XDRFileStream.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief XDRFileStream.cpp
 * @author Kai Buschulte
 * @date 15.04.2010
 * @since 1.0.0
 */

// hpp
#include <lama/io/XDRFileStream.hpp>

#include <sstream>

namespace lama
{

XDRFileStream::XDRFileStream( const std::string & filename, const std::ios_base::openmode openmode )
    : m_openmode( openmode )
{
    openFilePointer( filename, openmode );
}

XDRFileStream::~XDRFileStream()
{
    close();
}

bool XDRFileStream::is_open()
{
    return m_filestream.is_open();
}

void XDRFileStream::openFilePointer( const std::string filename, const std::ios_base::openmode openmode )
{

    if ( openmode & std::ios::in )
    {
        if ( openmode & std::ios::out )
        {
            throw Exception( "XDRFileStream: Stream can not be in in and out mode!" );
        }
    }
    else if ( openmode & std::ios::out )
    {
    }
    else
    {
        throw Exception( "Unknown file access mode for opening file pointer to access XDR streams." );
    }

    m_openmode = openmode | std::ios::binary;
    m_filename = filename;

    m_filestream.open( m_filename.c_str(), m_openmode );

    if ( !( m_filestream.is_open() ) )
    {
        throw Exception( ( std::string( "XDRFileStream: Can not open file!" ) + m_filename ).c_str() );
    }

}

void XDRFileStream::close()
{
    if ( is_open() )
    {
        m_filestream.close();
    }
}

bool XDRFileStream::isLittleEndian()
{
    int a = 1;
    char *ch = reinterpret_cast<char*>( &a );
    return static_cast<int>( *ch ) != 0;
}

} //namespace lama
