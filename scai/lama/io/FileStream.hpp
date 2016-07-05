/**
 * @file FileStream.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Expansion of the std::fstream class
 * @author Jan Ecker
 * @date 16.03.2016
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library

// scai internal libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/TemplateSpecifier.hpp>
#include <scai/common/macros/loop.hpp>

#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/hmemo.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/ReductionOp.hpp>


#include <fstream>
#include <string>
#include <typeinfo>
// TODO: remove
#include <iostream>

static inline void endian_convert( char* out, const char* in, const int n, const int size )
{
    std::cout << "endian_convert, n = " << n << ", size = " << size << std::endl;

    size_t offs = 0;

    for ( int i = 0; i < n; ++i )
    {
        char* out_s = out + offs;
        const char* in_s = in + offs;
        if ( size == 1 )
        {
            out_s[0] = in_s[0];
        }
        else if ( size == 2 )
        {
            out_s[0] = in_s[1];
            out_s[1] = in_s[0];
        }
        if ( size == 4 )
        {
            out_s[0] = in_s[3];
            out_s[1] = in_s[2];
            out_s[2] = in_s[1];
            out_s[3] = in_s[0];
        }
        else if ( size == 8 )
        {
            out_s[0] = in_s[7];
            out_s[1] = in_s[6];
            out_s[2] = in_s[5];
            out_s[3] = in_s[4];
            out_s[4] = in_s[3];
            out_s[5] = in_s[2];
            out_s[6] = in_s[1];
            out_s[7] = in_s[0];
        }
        else
        {
             COMMON_THROWEXCEPTION( "unsupported size" )
        }
        offs += size;
    }
}

// TODO: do initialization of all arrays using set() (correct first touch)
namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT FileStream : public std::fstream
{
public:

    typedef enum
    {
        BIG,
        LITTLE
    } Endian;

    FileStream( const std::string& filename, ios_base::openmode mode, Endian usedEndian = LITTLE );

    inline void open( const std::string& filename, ios_base::openmode mode, Endian usedEndian = LITTLE );

    template<typename ValueType>
    inline void write( const hmemo::HArray<ValueType>& data,
                       const ValueType offset,
                       const common::scalar::ScalarType type,
                       const char delimiter = ' ' );

    template<typename ValueType>
    inline void read(  hmemo::HArray<ValueType>& data,
                       const IndexType size,
                       const ValueType offset,
                       const common::scalar::ScalarType type,
                       const char delimiter = ' ' );

private:

    template<typename FileType, typename DataType>
    inline void _write( const hmemo::HArray<DataType>& data,
                        const FileType offset,
                        const char delimiter = ' ' );

    template<typename FileType, typename DataType>
    inline void _read(  hmemo::HArray<DataType>& data,
                        const IndexType size,
                        const DataType offset,
                        const char delimiter = ' ' );
    typedef enum
    {
        FORMATTED,
        BINARY
    } Type;

    Type mType;
    Endian mUsedEndian;

    /** Logger for this class */
    SCAI_LOG_DECL_STATIC_LOGGER( logger )


    static Endian mMachineEndian;
    static Endian _determineMachineEndian();

    template<typename ValueType, typename TList> struct Wrapper;

}; // class FileStream

template<typename ValueType>
struct FileStream::Wrapper<ValueType, common::mepr::NullType>
{
    static bool __write( FileStream&, const hmemo::HArray<ValueType>&, const ValueType, const common::scalar::ScalarType, const char )
    {
        return false;
    }

    static bool __read( FileStream&, hmemo::HArray<ValueType>&, const IndexType, const ValueType, const common::scalar::ScalarType, const char )
    {
        return false;
    }
};

template<typename ValueType, typename H, typename T>
struct FileStream::Wrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static bool __write( FileStream& fs, const hmemo::HArray<ValueType>& data, const ValueType offset, const common::scalar::ScalarType type, const char delimiter )
    {
        if ( type == common::TypeTraits<H>::stype )
        {
            fs._write<H, ValueType>( data, offset, delimiter );
            return true;
        }
        else
        {
            return Wrapper<ValueType, T>::__write( fs, data, offset, type, delimiter );
        }
    }

    static bool __read( FileStream& fs, hmemo::HArray<ValueType>& data, const IndexType size, const ValueType offset, const common::scalar::ScalarType type, const char delimiter )
    {
        if ( type == common::TypeTraits<H>::stype )
        {
            fs._read<H, ValueType>( data, size, offset, delimiter );
            return true;
        }
        else
        {
            return Wrapper<ValueType, T>::__read( fs, data, size, offset, type, delimiter );
        }
    }
};

inline FileStream::FileStream( const std::string& filename, ios_base::openmode mode, Endian usedEndian /* = LITTLE */ )
{
    open( filename, mode, usedEndian );
}

inline void FileStream::open( const std::string& filename, ios_base::openmode mode, Endian usedEndian /* = LITTLE */ )
{
    std::fstream::open( filename.c_str(), mode );

    if ( mode & std::ios::binary )
    {
        mType = BINARY;
    }
    else
    {
        mType = FORMATTED;
    }

    this->mUsedEndian = usedEndian;

    if ( !is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open file " + filename + "." )
    }
}

template<typename ValueType>
inline void FileStream::write( const hmemo::HArray<ValueType>& data,
                               const ValueType offset,
                               const common::scalar::ScalarType type,
                               const char delimiter /* = ' ' */ )
{
    if ( !Wrapper<ValueType, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::__write( *this, data, offset, type, delimiter ) )
    {
        switch ( type )
        {
            case common::scalar::INT:
                _write<int, ValueType>( data, offset, delimiter );
                break;

            case common::scalar::LONG:
                _write<long, ValueType>( data, offset, delimiter );
                break;

            case common::scalar::INDEX_TYPE:
                _write<IndexType, ValueType>( data, offset, delimiter );
                break;

            case common::scalar::INTERNAL:
                _write<ValueType, ValueType>( data, offset, delimiter );
                break;

            default:
                SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << type )
                break;
        }
    }
}

template<typename ValueType>
inline void FileStream::read( hmemo::HArray<ValueType>& data,
                              const IndexType size,
                              const ValueType offset,
                              const common::scalar::ScalarType type,
                              const char delimiter /* = ' ' */ )
{
    if ( !Wrapper<ValueType, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::__read( *this, data, size, offset, type, delimiter ) )
    {
        switch ( type )
        {
            case common::scalar::INT:
                _read<int, ValueType>( data, size, offset, delimiter );
                break;

            case common::scalar::LONG:
                _read<long, ValueType>( data, size, offset, delimiter );
                break;

            case common::scalar::INDEX_TYPE:
                _read<IndexType, ValueType>( data, size, offset, delimiter );
                break;

            case common::scalar::INTERNAL:
                _read<ValueType, ValueType>( data, size, offset, delimiter );
                break;

            default:
                SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << type )
                break;
        }
    }
}

template<typename FileType, typename DataType>
inline void FileStream::_write( const hmemo::HArray<DataType>& data,
                                const FileType offset,
                                const char delimiter /* = ' ' */ )
{
    SCAI_LOG_INFO( logger, "write array data <" << common::TypeTraits<DataType>::id() << "> to <"
                   << common::TypeTraits<FileType>::id() << ">, offset = " << offset )
    SCAI_ASSERT( is_open(), "FileStream is not opened" );

    if ( offset == 0 && typeid( FileType ) == typeid( DataType ) )
    {
        hmemo::ReadAccess<DataType> dataRead( data );

        switch ( mType )
        {
            case BINARY:
            {
                if ( mUsedEndian == mMachineEndian )
                {
                    std::fstream::write( reinterpret_cast<const char*>( dataRead.get() ), sizeof( DataType ) * data.size() );
                }
                else
                {
                    scai::common::scoped_array<DataType> tmp( new DataType[ data.size() ] );

                    endian_convert( reinterpret_cast<char*>( tmp.get() ), 
                                    reinterpret_cast<const char*>( dataRead.get() ), 
                                    data.size(), sizeof( DataType ) );

                    std::fstream::write( reinterpret_cast<const char*>( tmp.get() ), sizeof( DataType ) * data.size() );
                }

                break;
            }

            case FORMATTED:
            {
                const DataType* dataPtr = dataRead.get();

                for ( IndexType i = 0; i < data.size(); ++i )
                {
                    *this << dataPtr[i] << delimiter;
                }

                break;
            }
        }
    }
    else
    {
        hmemo::HArray<FileType> buffer;
        {
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<FileType, DataType> > set;
            static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<FileType> > setVal;
            hmemo::ContextPtr loc = data.getValidContext();
            setVal.getSupportedContext( loc );
            set.getSupportedContext( loc );
            hmemo::ReadAccess<DataType> dataRead( data, loc );
            hmemo::WriteOnlyAccess<FileType> bufferWrite( buffer, loc, data.size() );
            set[loc]( bufferWrite, dataRead, data.size(), utilskernel::reduction::COPY );
            setVal[loc]( bufferWrite, buffer.size(), offset, utilskernel::reduction::ADD );
        }
        hmemo::ReadAccess<FileType> bufferRead( buffer );

        switch ( mType )
        {
            case BINARY:
            {
                std::fstream::write( reinterpret_cast<const char*>( bufferRead.get() ), sizeof( FileType )*data.size() );
                break;
            }

            case FORMATTED:
            {
                const FileType* bufferPtr = bufferRead.get();

                for ( IndexType i = 0; i < data.size(); ++i )
                {
                    *this << bufferPtr[i] << delimiter;
                }

                break;
            }
        }
    }

    std::fstream::flush();
}

template<typename FileType, typename DataType>
inline void FileStream::_read( hmemo::HArray<DataType>& data,
                               const IndexType size,
                               const DataType offset,
                               const char delimiter /* = ' ' */ )
{
    SCAI_LOG_INFO( logger, "read array data <" << common::TypeTraits<FileType>::id() << "> to <"
                   << common::TypeTraits<DataType>::id() << ">, offset = " << offset << ", size = " << size )
    hmemo::WriteOnlyAccess<DataType> dataWrite( data, size );

    if ( typeid( FileType ) == typeid( DataType ) )
    {
        switch ( mType )
        {
            case BINARY:
            {
                if ( mUsedEndian == mMachineEndian )
                {
                    std::fstream::read( reinterpret_cast<char*>( dataWrite.get() ), sizeof( DataType ) * size );
                }
                else
                {
                    scai::common::scoped_array<DataType> tmp( new DataType[ data.size() ] );

                    std::fstream::read( reinterpret_cast<char*>( tmp.get() ), sizeof( DataType ) * size );

                    endian_convert( reinterpret_cast<char*>( dataWrite.get() ), 
                                    reinterpret_cast<const char*>( tmp.get() ), 
                                    size, sizeof( DataType ) );

                }


                if ( !*this )
                {
                    COMMON_THROWEXCEPTION( "Error reading data!" )
                }

                break;
            }

            case FORMATTED:
            {
                std::string buffer;
                std::stringstream ssBuffer;

                for ( int i = 0; i < size; ++i )
                {
                    std::getline( *this, buffer, delimiter );

                    if ( !*this )
                    {
                        COMMON_THROWEXCEPTION( "Unexpected end of file" )
                    }

                    ssBuffer.clear();
                    ssBuffer << buffer;
                    ssBuffer >> dataWrite[i];
                }

                break;
            }
        }
    }
    else
    {
        hmemo::HArray<FileType> buffer;
        hmemo::WriteOnlyAccess<FileType> bufferWrite( buffer, size );

        switch ( mType )
        {
            case BINARY:
            {
                std::fstream::read( reinterpret_cast<char*>( bufferWrite.get() ), sizeof( FileType )*size );

                if ( !*this )
                {
                    COMMON_THROWEXCEPTION( "Error reading data!" )
                }

                break;
            }

            case FORMATTED:
            {
                std::string buffer;
                std::stringstream ssBuffer;

                for ( int i = 0; i < size; ++i )
                {
                    std::getline( *this, buffer, delimiter );

                    if ( !*this )
                    {
                        COMMON_THROWEXCEPTION( "Unexpected end of file" )
                    }

                    ssBuffer.clear();
                    ssBuffer << buffer;
                    ssBuffer >> bufferWrite[i];
                }

                break;
            }
        }

        static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<DataType, FileType> > set;
        hmemo::ContextPtr loc = data.getValidContext();
        set.getSupportedContext( loc );
        set[loc]( dataWrite, bufferWrite, size, utilskernel::reduction::COPY );
    }

    if ( offset != 0 )
    {
        static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<DataType> > setVal;
        hmemo::ContextPtr loc = data.getValidContext();
        setVal.getSupportedContext( loc );
        setVal[loc]( dataWrite, size, offset, utilskernel::reduction::ADD );
    }
}

inline FileStream::Endian FileStream::_determineMachineEndian()
{
    int a = 1;
    char* ch = reinterpret_cast<char*>( &a );

    if ( static_cast<int>( *ch ) != 0 )
    {
        return LITTLE;
    }
    else
    {
        return BIG;
    }
}

} /* end namespace lama */

} /* end namespace scai */
