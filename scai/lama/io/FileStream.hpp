/**
 * @file FileStream.hpp
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
 * @brief Expansion of the std::fstream class
 * @author Jan Ecker
 * @date 16.03.2016
 * @since 2.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/TemplateSpecifier.hpp>
#include <scai/common/macros/typeloop.hpp>

#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/hmemo.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>


#include <fstream>
#include <string>
#include <typeinfo>
// TODO: remove
#include <iostream>


// TODO: do initialization of all arrays using set() (correct first touch)
namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT FileStream : public std::fstream
{
public:
    FileStream( const std::string& filename, ios_base::openmode mode, bool useLittleEndian = true );

    template<typename ValueType>
    void write( const hmemo::HArray<ValueType>& data,
                const int offset,
                const common::scalar::ScalarType type,
                const char delimiter = ' ' );

    template<typename ValueType>
    void read(  hmemo::HArray<ValueType>& data,
                const IndexType size,
                const int offset,
                const common::scalar::ScalarType type,
                const char delimiter = ' ' );

private:
    bool isBinary;
    bool isLittleEndian;
    bool useLittleEndian;

    /** Logger for this class */
    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename FileType, typename DataType>
    void _write( const hmemo::HArray<DataType>& data,
                 const int offset,
                 const char delimiter = ' ' );

    template<typename FileType, typename DataType>
    void _read(  hmemo::HArray<DataType>& data,
                 const IndexType size,
                 const int offset,
                 const char delimiter = ' ' );

    bool _isLittleEndian();

}; // class FileStream

SCAI_LOG_DEF_LOGGER( FileStream::logger, "IO.FileStream" )

FileStream::FileStream( const std::string& filename, ios_base::openmode mode, bool useLittleEndian /* = true */ )
{
    std::fstream::open( filename.c_str(), mode );
    if( mode & std::ios::binary )
    {
        isBinary = true;
    }
    else
    {
        isBinary = false;
    }
    this->useLittleEndian = useLittleEndian;
    isLittleEndian = _isLittleEndian();

    if( !is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open file " + filename + "." )
    }
}

template<typename ValueType>
void FileStream::write( const hmemo::HArray<ValueType>& data,
                        const int offset,
                        const common::scalar::ScalarType type,
                        const char delimiter /* = ' ' */ )
{
    switch(type){
    // special cases for handling IndexType, int and long, as these are not properly supported yet
    case common::scalar::INDEX_TYPE:
        _write<IndexType,ValueType>( data, offset, delimiter );
        break;
    case common::scalar::INT:
        _write<int,ValueType>( data, offset, delimiter );
        break;
    case common::scalar::LONG:
        _write<long,ValueType>( data, offset, delimiter );
        break;

    // generate cases for all scalar types
#define SCAI_LAMA_FILESTREAM_WRITE( _type )                \
    case ( common::TypeTraits<_type>::stype ):             \
        _write<_type,ValueType>( data, offset, delimiter ); \
        break;
    SCAI_COMMON_TYPELOOP( ARITHMETIC_HOST_CNT, SCAI_LAMA_FILESTREAM_WRITE, ARITHMETIC_HOST )
#undef SCAI_LAMA_FILESTREAM_WRITE

    case common::scalar::INTERNAL:
        _write<ValueType,ValueType>( data, offset, delimiter );
        break;
    default:
        SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << type )
        break;
    }
}

template<typename ValueType>
void FileStream::read(  hmemo::HArray<ValueType>& data,
                        const IndexType size,
                        const int offset,
                        const common::scalar::ScalarType type,
                        const char delimiter /* = ' ' */ )
{
    switch(type){
    // special cases for handling IndexType, int and long, as these are not properly supported yet
    case common::scalar::INDEX_TYPE:
        _read<IndexType,ValueType>( data, size, offset, delimiter );
        break;
    case common::scalar::INT:
        _read<int,ValueType>( data, size, offset, delimiter );
        break;
    case common::scalar::LONG:
        _read<long,ValueType>( data, size, offset, delimiter );
        break;

    // generate cases for all scalar types
#define SCAI_LAMA_FILESTREAM_READ( _type )                       \
    case ( common::TypeTraits<_type>::stype ):                   \
        _read<_type,ValueType>( data, size, offset, delimiter ); \
        break;
    SCAI_COMMON_TYPELOOP( ARITHMETIC_HOST_CNT, SCAI_LAMA_FILESTREAM_READ, ARITHMETIC_HOST )
#undef SCAI_LAMA_FILESTREAM_READ

    case common::scalar::INTERNAL:
        _read<ValueType,ValueType>( data, size, offset, delimiter );
        break;
    default:
        SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << type )
        break;
    }
}

template<typename FileType, typename DataType>
void FileStream::_write( const hmemo::HArray<DataType>& data,
                         const int offset,
                         const char delimiter /* = ' ' */ )
{
    SCAI_LOG_INFO( logger, "write array data <" << common::TypeTraits<DataType>::id() << "> to <"
                           << common::TypeTraits<FileType>::id() << ">, offset = " << offset )

    assert( is_open() );


    if( offset == 0 && typeid(FileType) == typeid(DataType) )
    {
        hmemo::ReadAccess<DataType> dataRead( data );
        if ( isBinary )
        {
            if( useLittleEndian == isLittleEndian )
            {
                std::fstream::write( reinterpret_cast<const char*>( dataRead.get() ), sizeof(DataType)*data.size() );
            }
            else
            {
                COMMON_THROWEXCEPTION( "Error, wrong Endian!!!" )
            }
        }
        else
        {
            const DataType* dataPtr = dataRead.get();
            for( IndexType i = 0; i < data.size(); ++i )
            {
                *this << dataPtr[i] << delimiter;
            }
        }
    }
    else
    {
        hmemo::HArray<FileType> buffer;
        {
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<FileType, DataType> > set;
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::addScalar<FileType> > addScalar;
            hmemo::ContextPtr loc = addScalar.getValidContext(set.getValidContext(data.getValidContext()));

            hmemo::ReadAccess<DataType> dataRead( data, loc );
            hmemo::WriteOnlyAccess<FileType> bufferWrite( buffer, loc, data.size() );

            set[loc](bufferWrite, dataRead, data.size(), common::reduction::COPY);
            addScalar[loc](bufferWrite, buffer.size(), offset);
        }
        hmemo::ReadAccess<FileType> bufferRead( buffer );

        if( isBinary )
        {
            std::fstream::write( reinterpret_cast<const char*>( bufferRead.get() ), sizeof(FileType)*data.size() );
        }
        else
        {
            const FileType* bufferPtr = bufferRead.get();
            for( IndexType i = 0; i < data.size(); ++i )
            {
                *this << bufferPtr[i] << delimiter;
            }
        }
    }
    std::fstream::flush();
}

template<typename FileType, typename DataType>
void FileStream::_read( hmemo::HArray<DataType>& data,
                        const IndexType size,
                        const int offset,
                        const char delimiter /* = ' ' */ )
{
    SCAI_LOG_INFO( logger, "read array data <" << common::TypeTraits<FileType>::id() << "> to <"
                           << common::TypeTraits<DataType>::id() << ">, offset = " << offset << ", size = " << size )

    hmemo::WriteOnlyAccess<DataType> dataWrite( data, size );

    if(typeid(FileType) == typeid(DataType) )
    {
        if( isBinary )
        {
            std::fstream::read( reinterpret_cast<char*>( dataWrite.get() ), sizeof(DataType)*size );
            if ( !*this )
            {
                COMMON_THROWEXCEPTION( "Error reading data!" )
            }
        }
        else
        {
            std::string buffer;
            std::stringstream ssBuffer;
            for( int i=0; i < size; ++i )
            {
                std::getline(*this, buffer, delimiter );
                if ( !*this )
                {
                    COMMON_THROWEXCEPTION( "Unexpected end of file" )
                }

                ssBuffer.clear();
                ssBuffer << buffer;
                ssBuffer >> dataWrite[i];
            }
        }
    }
    else
    {
        hmemo::HArray<FileType> buffer;
        hmemo::WriteOnlyAccess<FileType> bufferWrite( buffer, size );

        if( isBinary )
        {
            std::fstream::read( reinterpret_cast<char*>( bufferWrite.get() ), sizeof(FileType)*size );
            if ( !*this )
            {
                COMMON_THROWEXCEPTION( "Error reading data!" )
            }
        }
        else
        {
            std::string buffer;
            std::stringstream ssBuffer;
            for( int i=0; i < size; ++i )
            {
                std::getline(*this, buffer, delimiter );
                if ( !*this )
                {
                    COMMON_THROWEXCEPTION( "Unexpected end of file" )
                }

                ssBuffer.clear();
                ssBuffer << buffer;
                ssBuffer >> bufferWrite[i];
            }
        }


        static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<DataType, FileType> > set;
        hmemo::ContextPtr loc = set.getValidContext(data.getValidContext());
        set[loc](dataWrite, bufferWrite, size, common::reduction::COPY);
    }

    if( offset != 0 ){
        static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::addScalar<DataType> > addScalar;
        hmemo::ContextPtr loc = addScalar.getValidContext(data.getValidContext());
        addScalar[loc](dataWrite, size, offset);
    }
}

bool FileStream::_isLittleEndian()
{
    int a = 1;
    char *ch = reinterpret_cast<char*>( &a );
    return static_cast<int>( *ch ) != 0;
}


} /* end namespace lama */

} /* end namespace scai */
