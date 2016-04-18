/*
 * XDRUtils.h
 *
 *  Created on: Mar 8, 2016
 *      Author: eschricker
 */

#pragma once

// local library

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

    /** Reading binary data with implicit conversion from FileType to DataType and adding offset */
    template<typename FileDataType,typename UserDataType>
    static void readBinaryData( std::fstream& inFile, UserDataType data[], const IndexType n, const IndexType offset = 0 );

    /** Writing binary data with implicit conversion from DataType to FileType and adding offset */
    template<typename FileType, typename DataType>
    static void writeBinary( std::fstream& outFile, const DataType data[], const IndexType n, const IndexType offset = 0 );
};

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
