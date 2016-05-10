/**
 * @file lama/io/IOUtils.hpp
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
 * @brief ToDo: Missing description in ./lama/io/IOUtils.hpp
 * @author eschricker
 * @date 16.03.2015
 */

#pragma once

// local library
#include <scai/lama/io/XDRFileStream.hpp>

// scai internal libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/unique_ptr.hpp>

// others
#include <ostream>
#include <fstream>
#include <typeinfo>

namespace scai
{

namespace lama
{

class IOUtils
{
public:

    /** Own data type for file size, int might not be sufficient */
    typedef size_t file_size_t;

    /** Get the size of a file */
    static file_size_t getFileSize( const char* filename );

    template<typename FileType, typename DataType>
    static void readXDR( XDRFileStream& inFile, DataType* data, const IndexType n, const IndexType offset = 0 );

    template<typename FileType, typename DataType>
    static void writeXDR( XDRFileStream& outFile, const DataType* data, const IndexType n, const IndexType offset = 0 );

    /** Reading binary data with implicit conversion from FileType to DataType and adding offset */
    template<typename FileDataType,typename UserDataType>
    static void readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n, const IndexType offset = 0 );

    /** Writing binary data with implicit conversion from DataType to FileType and adding offset */
    template<typename FileType, typename DataType>
    static void writeBinary( std::fstream& outFile, const DataType data[], const IndexType n, const IndexType offset = 0 );
};

template<typename FileType, typename DataType>
void IOUtils::readXDR( XDRFileStream& inFile, DataType* data, const IndexType n, const IndexType offset )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        inFile.read( data, n );
        return;
    }

    // allocate a temporary buffer for n values of FileDataType to read the data

    common::scoped_array<FileType> buffer( new FileType[n] );

    inFile.read( buffer.get(), n );

    for( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<DataType>( buffer[i] + offset );
    }
}

template<typename FileType, typename DataType>
void IOUtils::writeXDR( XDRFileStream& outFile, const DataType* data, const IndexType n, const IndexType offset )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        outFile.write( data, n ); // no conversion needed
        return;
    }

    // so user data has to be converted in file type data

    common::scoped_array<FileType> buffer( new FileType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] + offset );
    }

    outFile.write( buffer.get(), n );
}

template<typename FileType,typename DataType>
void IOUtils::readBinaryData( std::fstream& inFile, DataType data[], const IndexType n, const IndexType offset )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        inFile.read( reinterpret_cast<char*>( data ), sizeof(DataType) * n );
    }
    else
    {
        // allocate buffer with file data for type conversion

        common::scoped_array<FileType> buffer( new FileType[n] );

        inFile.read( reinterpret_cast<char*>( buffer.get() ), sizeof(FileType) * n );

        for( IndexType i = 0; i < n; i++ )
        {
            data[i] = static_cast<DataType>( buffer[i] + offset );
        }
    }
}

template<typename FileType,typename DataType>
void IOUtils::writeBinary( std::fstream& outFile, const DataType data[], const IndexType n, const IndexType offset )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        outFile.write( reinterpret_cast<const char*>( data ), sizeof(DataType) * n );
        outFile.flush();
        return;
    }

    // allocate buffer for type conversion

    common::scoped_array<FileType> buffer( new FileType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] + offset );
    }

    outFile.write( reinterpret_cast<const char*>( buffer.get() ), sizeof(FileType) * n );
    outFile.flush();
    return;
}

} /* end namespace lama */

} /* end namespace scai */
