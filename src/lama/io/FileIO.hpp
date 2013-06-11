/**
 * @file FileIO.hpp
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
 * @brief Class providing IO routines for matrix storage
 * @author Thomas Brandes
 * @date 31.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_FILE_IO_HPP_
#define LAMA_FILE_IO_HPP_

// for dll_import
#include <lama/config.hpp>

#include <lama/LAMATypes.hpp>

#include <boost/scoped_array.hpp>
#include <fstream>
#include <typeinfo>

namespace lama
{

/** This class provides File operations needed for LAMA.  */

class LAMA_DLL_IMPORTEXPORT FileIO
{

public:

    /** Own data type for file size, int might not be sufficient */

    typedef size_t file_size_t;

    /** Get the size of a file */

    static file_size_t getFileSize( const char* filename );

    /** Reading binary data with implicit conversion from FileDataType to UserDataType */

    template<typename FileDataType, typename UserDataType>
    static void readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n );
  
    /** Reading binary data with implicit conversion from FileDataType to UserDataType and adding offset */

    template<typename FileDataType, typename UserDataType, int offset>
    static void readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n );
};

/* -------------------------------------------------------------------------- */

template<typename FileDataType, typename UserDataType>
void FileIO::readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n )
{
    if ( typeid( FileDataType ) == typeid( UserDataType ) )
    {
        // no type conversion needed

        inFile.read( reinterpret_cast<char*>( data ), sizeof( UserDataType ) * n );
        return;
    }

    // allocate buffer with file data for type conversion

    boost::scoped_array<FileDataType> buffer( new FileDataType[n] );

    inFile.read( reinterpret_cast<char*>( buffer.get() ), sizeof( FileDataType ) * n );

    for ( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<UserDataType>( buffer[i] );
    }
}

template<typename FileDataType, typename UserDataType, int offset>
void FileIO::readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n )
{
    if ( offset == 0 )
    {
        readBinaryData<FileDataType, UserDataType>( inFile, data, n );
        return;
    }

    // allocate buffer for type conversion and/or adding offset

    boost::scoped_array<FileDataType> buffer( new FileDataType[n] );

    inFile.read( reinterpret_cast<char*>( buffer.get() ), sizeof( FileDataType ) * n );

    for ( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<UserDataType>( buffer[i] + offset );
    }
}

}
#endif // LAMA_FILE_IO_HPP_
