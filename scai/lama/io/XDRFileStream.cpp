/**
 * @file XDRFileStream.cpp
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
 * @brief XDRFileStream.cpp
 * @author Kai Buschulte
 * @date 15.04.2010
 */

// hpp
#include <scai/lama/io/XDRFileStream.hpp>

// std
#include <sstream>

namespace scai
{

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

    if( openmode & std::ios::in )
    {
        if( openmode & std::ios::out )
        {
            COMMON_THROWEXCEPTION( "XDRFileStream: Stream can not be in in and out mode!" );
        }
    }
    else if( openmode & std::ios::out )
    {
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unknown file access mode for opening file pointer to access XDR streams." );
    }

    m_openmode = openmode | std::ios::binary;
    m_filename = filename;

    m_filestream.open( m_filename.c_str(), m_openmode );

    if( !( m_filestream.is_open() ) )
    {
        COMMON_THROWEXCEPTION( "XDRFileStream: Can not open file " << m_filename );
    }

}

void XDRFileStream::close()
{
    if( is_open() )
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

} /* end namespace lama */

} /* end namespace scai */
