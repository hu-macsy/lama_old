/**
 * @file CollectiveFile.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Collectice/concurrent file I/O 
 * @author Thomas Brandes
 * @date 14.01.2019
 */

#include <scai/dmemo/CollectiveFile.hpp>

#include <scai/common/macros/loop.hpp>

namespace scai
{

namespace dmemo
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CollectiveFile::logger, "CollectiveFile" )

/* -------------------------------------------------------------------------- */

CollectiveFile::CollectiveFile( CommunicatorPtr comm ) :

    mComm( comm ),
    mOffset( 0 )
{
}

CollectiveFile::~CollectiveFile()
{
}

/* -------------------------------------------------------------------------- */

void CollectiveFile::set( const char* filename, size_t offset )
{
    mFileName = filename;
    mOffset = offset;
}

/* -------------------------------------------------------------------------- */
/*  I/O conversion wrapper class                                              */
/* -------------------------------------------------------------------------- */

template<typename ValueType, typename TList>
struct IOWrapper;

template<typename ValueType>
struct IOWrapper<ValueType, common::mepr::NullType>
{
    static void writeSingle( CollectiveFile&, const ValueType[], const IndexType, const common::ScalarType fileType )
    {
        COMMON_THROWEXCEPTION( "Unsupported file type: " << fileType );
    }

    static void readSingle( CollectiveFile&, ValueType[], const IndexType, const common::ScalarType fileType )
    {
        COMMON_THROWEXCEPTION( "Unsupported file type: " << fileType );
    }

    static void writeAll( CollectiveFile&, const ValueType[], const IndexType, const IndexType, const common::ScalarType fileType )
    {
        COMMON_THROWEXCEPTION( "Unsupported file type: " << fileType );
    }

    static void readAll( CollectiveFile&, ValueType[], const IndexType, const IndexType, const common::ScalarType fileType )
    {
        COMMON_THROWEXCEPTION( "Unsupported file type: " << fileType );
    }
};

template<typename ValueType, typename FileType, typename TailTypes>
struct IOWrapper<ValueType, common::mepr::TypeList<FileType, TailTypes> >
{
    static void writeSingle( CollectiveFile& file, const ValueType array[], const IndexType n, const common::ScalarType fileType )
    {
        if ( common::TypeTraits<FileType>::stype == fileType )
        {
            std::unique_ptr<FileType[]> tmpArray( new FileType[n] );

            for ( IndexType i = 0; i < n; ++i )
            {
                tmpArray[i] = static_cast<FileType>( array[i] );
            }

            file.writeSingle( tmpArray.get(), n );
        }
        else
        {
            IOWrapper<ValueType, TailTypes>::writeSingle( file, array, n, fileType );
        }
    }

    static void readSingle( CollectiveFile& file, ValueType array[], const IndexType n, const common::ScalarType fileType )
    {
        if ( common::TypeTraits<FileType>::stype == fileType )
        {
            std::unique_ptr<FileType[]> tmpArray( new FileType[n] );

            file.readSingle( tmpArray.get(), n );

            for ( IndexType i = 0; i < n; ++i )
            {
                array[i] = static_cast<ValueType>( tmpArray[i] );
            }
        }
        else
        {
            IOWrapper<ValueType, TailTypes>::readSingle( file, array, n, fileType );
        }
    }

    static void writeAll( CollectiveFile& file, const ValueType array[], const IndexType n, const IndexType offset, const common::ScalarType fileType )
    {
        if ( common::TypeTraits<FileType>::stype == fileType )
        {
            std::unique_ptr<FileType[]> tmpArray( new FileType[n] );

            for ( IndexType i = 0; i < n; ++i )
            {
                tmpArray[i] = static_cast<FileType>( array[i] );
            }

            file.writeAll( tmpArray.get(), n, offset );
        }
        else
        {
            IOWrapper<ValueType, TailTypes>::writeAll( file, array, n, offset, fileType );
        }
    }

    static void readAll( CollectiveFile& file, ValueType array[], const IndexType n, const IndexType offset, const common::ScalarType fileType )
    {
        if ( common::TypeTraits<FileType>::stype == fileType )
        {
            std::unique_ptr<FileType[]> tmpArray( new FileType[n] );

            file.readAll( tmpArray.get(), n, offset );

            for ( IndexType i = 0; i < n; ++i )
            {
                array[i] = static_cast<ValueType>( tmpArray[i] );
            }
        }
        else
        {
            IOWrapper<ValueType, TailTypes>::readAll( file, array, n, offset, fileType );
        }
    }
};

/* -------------------------------------------------------------------------- */
/*  Read                                                                      */
/* -------------------------------------------------------------------------- */


template<typename ValueType>
void CollectiveFile::readSingle( ValueType array[], const IndexType n )
{
    const PartitionId MASTER = 0;

    if ( mComm->getRank() == MASTER )
    {
        auto stype = common::TypeTraits<ValueType>::stype;
        this->readSingleImpl( array, n, mOffset, stype );
    }

    mOffset += n * sizeof( ValueType );

    mComm->bcast( array, n, MASTER );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( ValueType array[], const IndexType n, const common::ScalarType fileType )
{
    if ( fileType == common::ScalarType::INTERNAL || fileType == common::TypeTraits<ValueType>::stype )
    {
        readSingle( array, n );
    }
    else
    {
        IOWrapper<ValueType, SCAI_ALL_TYPES_LIST>::readSingle( *this, array, n, fileType );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( ValueType& val, const common::ScalarType fileType )
{
    // consider val as a one-dimensional array with one entry

    this->readSingle( &val, 1, fileType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( hmemo::HArray<ValueType>& array, const IndexType n, const common::ScalarType fileType )
{
    auto wArray = hmemo::hostWriteOnlyAccess( array, n );
    this->readSingle( wArray.get(), n, fileType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readAll( ValueType local[], const IndexType n, const IndexType offset )
{
    // no more type conversion, read local array from local offset

    auto arrayType = common::TypeTraits<ValueType>::stype;

    this->readAllImpl( local, n, mOffset + offset * sizeof( ValueType ), arrayType );

    mOffset += mComm->sum( n ) * sizeof( ValueType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readAll( ValueType array[], const IndexType n, const IndexType offset, const common::ScalarType fileType )
{
    auto arrayType = common::TypeTraits<ValueType>::stype;

    if ( fileType == common::ScalarType::INTERNAL || fileType == arrayType )
    {
        readAll( array, n, offset );
    }
    else
    {
        IOWrapper<ValueType, SCAI_ALL_TYPES_LIST>::readAll( *this, array, n, offset, fileType );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readAll( hmemo::HArray<ValueType>& local, const IndexType size, const IndexType offset, const common::ScalarType fileType )
{
    auto wLocal = hmemo::hostWriteOnlyAccess( local, size );
    readAll( wLocal.get(), size, offset, fileType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readAll( hmemo::HArray<ValueType>& local, const IndexType size, const common::ScalarType fileType )
{
    IndexType offset = mComm->scan( size ) - size;
    readAll( local, size, offset, fileType );
}

/* -------------------------------------------------------------------------- */
/*   writeSingle                                                              */
/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const ValueType array[], const IndexType n )
{
    if ( mComm->getRank() == 0 )
    {
        auto stype = common::TypeTraits<ValueType>::stype;
        writeSingleImpl( mOffset, array, static_cast<size_t>( n ), stype );
    }

    mOffset += sizeof( ValueType ) * n;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const ValueType array[], const IndexType n, const common::ScalarType fileType )
{
    auto arrayType = common::TypeTraits<ValueType>::stype;

    if ( fileType == common::ScalarType::INTERNAL || fileType == arrayType )
    {
        writeSingle( array, n );
    }
    else
    {
        IOWrapper<ValueType, SCAI_ALL_TYPES_LIST>::writeSingle( *this, array, n, fileType );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const ValueType val, const common::ScalarType fileType )
{
    writeSingle( &val, 1, fileType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const hmemo::HArray<ValueType>& array, const common::ScalarType fileType )
{
    auto rLocal = hmemo::hostReadAccess( array );
    writeSingle( rLocal.get(), rLocal.size(), fileType );
}

/* -------------------------------------------------------------------------- */
/*   writeAll                                                                 */
/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const ValueType local[], const IndexType n, const IndexType offset )
{
    auto arrayType = common::TypeTraits<ValueType>::stype;

    writeAllImpl( mOffset + offset * sizeof( ValueType ), local, static_cast<size_t>( n ), arrayType );

    mOffset += mComm->sum( n ) * sizeof( ValueType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const ValueType array[], const IndexType n, const IndexType offset, const common::ScalarType fileType )
{
    auto arrayType = common::TypeTraits<ValueType>::stype;

    if ( fileType == common::ScalarType::INTERNAL || fileType == arrayType )
    {
        writeAll( array, n, offset );
    }
    else
    {
        IOWrapper<ValueType, SCAI_ALL_TYPES_LIST>::writeAll( *this, array, n, offset, fileType );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const hmemo::HArray<ValueType>& local, IndexType offset, const common::ScalarType fileType )
{
    auto rLocal = hmemo::hostReadAccess( local );
    writeAll( rLocal.get(), local.size(), offset, fileType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const hmemo::HArray<ValueType>& local, const common::ScalarType fileType )
{
    // scan communication needed to get for each processor its offset

    IndexType size = local.size();
    IndexType offset = mComm->scan( size );   // incluse scan to get offset

    this->writeAll( local, offset - size, fileType );
}

/* -------------------------------------------------------------------------- */

#define SCAI_COLLECTIVE_FILE_METHOD_INSTANTIATIONS( _type )       \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readSingle(                              \
        _type array[],                                            \
        const IndexType );                                        \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readSingle(                              \
        _type array[],                                            \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readSingle(                              \
        hmemo::HArray<_type>&,                                    \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readSingle(                              \
        _type& val,                                               \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readAll(                                 \
        _type array[],                                            \
        const IndexType,                                          \
        const IndexType );                                        \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readAll(                                 \
        _type array[],                                            \
        const IndexType,                                          \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readAll(                                 \
        hmemo::HArray<_type>&,                                    \
        const IndexType,                                          \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::readAll(                                 \
        hmemo::HArray<_type>&,                                    \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeSingle(                             \
        const _type array[],                                      \
        const IndexType );                                        \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeSingle(                             \
        const _type array[],                                      \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeSingle(                             \
        const _type val,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeSingle(                             \
        const hmemo::HArray<_type>& array,                        \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeAll(                                \
        const _type array[],                                      \
        const IndexType,                                          \
        const IndexType );                                        \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeAll(                                \
        const _type array[],                                      \
        const IndexType,                                          \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeAll(                                \
        const hmemo::HArray<_type>& array,                        \
        const IndexType,                                          \
        const common::ScalarType );                               \
                                                                  \
    template COMMON_DLL_IMPORTEXPORT                              \
    void CollectiveFile::writeAll(                                \
        const hmemo::HArray<_type>& array,                        \
        const common::ScalarType );                               \
                                                                  
SCAI_COMMON_LOOP( SCAI_COLLECTIVE_FILE_METHOD_INSTANTIATIONS, SCAI_ALL_TYPES )

}

}
