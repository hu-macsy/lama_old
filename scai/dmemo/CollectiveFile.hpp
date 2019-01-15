/**
 * @file CollectiveFile.hpp
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

#pragma once

#include <scai/logging.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/hmemo.hpp>

namespace scai
{

namespace dmemo
{

/**
 *  Abstract base class for a collective file, i.e. a file where all 
 *  processors of a communicator can write concurrently into it
 *  and read from it.
 */
class COMMON_DLL_IMPORTEXPORT CollectiveFile

{
public:

    ~CollectiveFile();

    /**
     *  @brief Pure method to open a file.
     *
     *  @param[in] fileName specifies the name of the file
     *  @param[in] fileMode is either "r" for reading a file or "w" for writing into it
     *
     *  This method might throw an IOException if the file could not be opened.
     */
    virtual void open( const char* fileName, const char* fileMode ) = 0;

    /**
     *  @brief Pure method to close the currently open file.
     */
    virtual void close() = 0;

    /**
     *  @brief Untyped version of writing data to the collective file by a single processor only
     *
     *  @param[in] offset specifes the fie position where to write
     *  @param[in] val    is pointer of the data to write
     *  @param[in] n      is the number of entries to write
     *  @param[in] stype  specifies the type of the data
     */
    virtual void writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype ) = 0;

    virtual void writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype ) = 0;

    /** 
     *  @brief One single processor writes an array of values at current file position
     *
     *  This method must be called by all processors even if only one single processor writes its values.
     */
    template<typename ValueType>
    void writeSingle( const ValueType val[], const IndexType n );

    /**
     *  @brief More convenient method for writing a single value only. 
     */
    template<typename ValueType>
    void writeSingle( const ValueType val );

    template<typename ValueType>
    void writeSingle( const hmemo::HArray<ValueType>& array );

    /** 
     *  @brief Each processor writes an array of values at current file position with an individual offset
     */
    template<typename ValueType>
    void writeAll( const hmemo::HArray<ValueType>& local, IndexType offset );

    /** 
     *  @brief Each processor writes an array of values at current file position with an individual offset
     */
    template<typename ValueType>
    void writeAll( const hmemo::HArray<ValueType>& local );

    virtual void readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype ) = 0;

    template<typename ValueType>
    void readSingle( ValueType val[], const IndexType n );

    template<typename ValueType>
    void readSingle( ValueType& val );

    template<typename ValueType>
    void readSingle( hmemo::HArray<ValueType>& val, const IndexType n );

    virtual void readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype ) = 0;

    template<typename ValueType>
    void readAll( hmemo::HArray<ValueType>& local, const IndexType size, const IndexType offset );

protected:

    void set( const char* filename, size_t offset );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    CollectiveFile( CommunicatorPtr comm );

    std::shared_ptr<const Communicator> mComm;

    std::string mFileName;

    size_t mOffset;    // current file position
};

/* -------------------------------------------------------------------------- */
/*   Implementation of template methods                                       */
/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const ValueType val[], const IndexType n )
{
    if ( mComm->getRank() == 0 )
    {
        auto stype = common::TypeTraits<ValueType>::stype;
        writeSingleImpl( mOffset, val, static_cast<size_t>( n ), stype );
    }

    mOffset += sizeof( ValueType ) * n;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const ValueType val )
{
    writeSingle( &val, 1 );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeSingle( const hmemo::HArray<ValueType>& array )
{
    auto rLocal = hmemo::hostReadAccess( array );
    writeSingle( rLocal.get(), rLocal.size() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const hmemo::HArray<ValueType>& local, IndexType offset )
{
    auto stype = common::TypeTraits<ValueType>::stype;
    auto rLocal = hmemo::hostReadAccess( local );
    this->writeAllImpl( mOffset + offset * sizeof( ValueType ), rLocal.get(), rLocal.size(), stype );

    mOffset += mComm->sum( local.size() ) * sizeof( ValueType );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::writeAll( const hmemo::HArray<ValueType>& local )
{
    // scan communication needed to get for each processor its offset

    IndexType size = local.size();
    IndexType offset = mComm->scan( size );   // incluse scan to get offset

    this->writeAll( local, offset - size );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( ValueType val[], const IndexType n )
{
    const PartitionId MASTER = 0;

    if ( mComm->getRank() == MASTER )
    {
        auto stype = common::TypeTraits<ValueType>::stype;
        this->readSingleImpl( val, n, mOffset, stype );
    }

    mOffset += n * sizeof( ValueType );

    mComm->bcast( val, n, MASTER );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( ValueType& val )
{
    this->readSingle( &val, 1 );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readSingle( hmemo::HArray<ValueType>& array, const IndexType n )
{
    auto wArray = hmemo::hostWriteOnlyAccess( array, n );
    this->readSingle( wArray.get(), n );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveFile::readAll( hmemo::HArray<ValueType>& local, const IndexType size, const IndexType offset )
{
    auto stype = common::TypeTraits<ValueType>::stype;
    auto wLocal = hmemo::hostWriteOnlyAccess( local, size );

    this->readAllImpl( wLocal.get(), size, mOffset + offset * sizeof( ValueType ), stype );

    mOffset += mComm->sum( local.size() ) * sizeof( ValueType );
}

/* -------------------------------------------------------------------------- */

}

}
