/*
 * XDRUtils.h
 *
 *  Created on: Mar 8, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/SCAITypes.hpp>
#include <scai/common/unique_ptr.hpp>

#include <ostream>

namespace scai
{

namespace lama
{

class IOUtils
{
public:

    template<typename FileType, typename DataType>
    static void readXDR( XDRFileStream& inFile, DataType* data, const IndexType n );

    template<typename FileType, typename DataType>
    static void writeXDR( XDRFileStream& outFile, const DataType* data, const IndexType n, const IndexType offset = 0 );

    template<typename FileType, typename DataType>
    static void writeBinary( std::fstream& outFile, const DataType data[], const IndexType n, const IndexType offset = 0 );
};

template<typename FileType, typename DataType>
void IOUtils::readXDR( XDRFileStream& inFile, DataType* data, const IndexType n )
{
    if( typeid(FileType) == typeid(DataType) )
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
        data[i] = static_cast<DataType>( buffer[i] );
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
