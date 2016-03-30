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
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/TemplateSpecifier.hpp>
#include <scai/common/macros/typeloop.hpp>

#include <scai/hmemo.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>


#include <fstream>
#include <string>
#include <typeinfo>
// TODO: remove
#include <iostream>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT FileStream : std::fstream
{
public:
    FileStream( const std::string& filename, ios_base::openmode mode = ios_base::in|ios_base::out );

    //void open( const std::string& filename );

    //void close();

    template<typename ValueType>
    void write( const hmemo::HArray<ValueType>& data, int offset, common::scalar::ScalarType );

protected:
    /** Logger for this class */
    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename FileType, typename DataType>
    void write( const hmemo::HArray<DataType>& data, int offset );

private:
    //std::fstream mFileStream;

    SCAI_DECLARE_TEMPLATESPECIFIER( SpecifierVO, template<typename ValueType, typename OtherValueType> )

    class Guard
    {
    public:
        Guard();
    };

    static Guard guard;

}; // class FileStream

SCAI_LOG_DEF_LOGGER( FileStream::logger, "FileStream" )

FileStream::FileStream( const std::string& filename, ios_base::openmode mode /* = ios_base::in|ios_base::out */ )
{
    std::fstream::open( filename.c_str(), mode );
}

//void FileStream::open( const std::string& filename )
//{
//    mFileStream.open( filename.c_str(), std::ios::out | std::ios::binary );
//
//}

//void FileStream::close()
//{
//    mFileStream.close();
//}

template<typename ValueType>
void FileStream::write( const hmemo::HArray<ValueType>& data, int offset, common::scalar::ScalarType type )
{
    switch(type){
    // special cases for handling IndexType, int and long, as these are not properly supported yet
    case common::scalar::INDEX_TYPE:
        write<IndexType,ValueType>( data, offset );
        break;
    case common::scalar::INT:
        write<int,ValueType>( data, offset );
        break;
    case common::scalar::LONG:
        write<long,ValueType>( data, offset );
        break;

    // generate cases for all scalar types
#define SCAI_LAMA_FILESTREAM_WRITE( _type )                \
    case ( common::TypeTraits<_type>::stype ):             \
        write<_type,ValueType>( data, offset );            \
        break;
    SCAI_COMMON_TYPELOOP( ARITHMETIC_HOST_CNT, SCAI_LAMA_FILESTREAM_WRITE, ARITHMETIC_HOST )
#undef SCAI_LAMA_FILESTREAM_WRITE

    case common::scalar::INTERNAL:
        write<ValueType,ValueType>( data, offset );
        break;
    default:
        SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << type )
        break;
    }
}

template<typename FileType, typename DataType>
void FileStream::write( const hmemo::HArray<DataType>& data, int offset )
{
    SCAI_LOG_INFO( logger, "write array data<" << common::TypeTraits<DataType>::id() << "> to <"
                           << common::TypeTraits<FileType>::id() << ">, offset = " << offset )

    assert( is_open() );

    if( offset == 0 && typeid(FileType) == typeid(DataType) ){
        hmemo::ReadAccess<DataType> dataRead( data );
        std::fstream::write( reinterpret_cast<const char*>( dataRead.get() ), sizeof(DataType)*data.size() );
        std::fstream::flush();
    }else{
        hmemo::HArray<FileType> buffer;
        {
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<FileType, DataType> > set;
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::addScalar<FileType> > addScalar;
            hmemo::ContextPtr loc = addScalar.getValidContext(set.getValidContext(data.getValidContext()));

            std::cout << "context: " << *loc << std::endl;

            hmemo::ReadAccess<DataType> dataRead( data, loc );
            hmemo::WriteOnlyAccess<FileType> bufferWrite( buffer, loc, data.size() );

            set[loc](bufferWrite, dataRead, data.size(), common::reduction::COPY);
            addScalar[loc](bufferWrite, buffer.size(), offset);
        }
        hmemo::ReadAccess<FileType> bufferRead( buffer );

        std::fstream::write( reinterpret_cast<const char*>( bufferRead.get() ), sizeof(FileType)*data.size() );
        std::fstream::flush();
    }
}


} /* end namespace lama */

} /* end namespace scai */
