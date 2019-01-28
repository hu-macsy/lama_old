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
 *  Abstract base class for a collective RAW file, i.e. a file where all 
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

    template<typename ValueType>
    void writeSingle( const ValueType array[], const IndexType n );

    /** 
     *  @brief One single processor writes an array of values at current file position
     *
     *  This method must be called by all processors even if only one single processor writes its values.
     */
    template<typename ValueType>
    void writeSingle( const ValueType array[], const IndexType n, const common::ScalarType fileType );

    /**
     *  @brief More convenient method for writing a single value only. 
     */
    template<typename ValueType>
    void writeSingle( const ValueType val, const common::ScalarType fileType = common::ScalarType::INTERNAL );

    template<typename ValueType>
    void writeSingle( const hmemo::HArray<ValueType>& array, const common::ScalarType fileType = common::ScalarType::INTERNAL );

    template<typename ValueType>
    void writeAll( const ValueType local[], const IndexType n, const IndexType offset );

    template<typename ValueType>
    void writeAll( const ValueType local[], const IndexType n, const IndexType offset, const common::ScalarType fileType );

    /** 
     *  @brief Each processor writes an array of values at current file position with an individual offset
     */
    template<typename ValueType>
    void writeAll( const hmemo::HArray<ValueType>& local, const IndexType offset, 
                   const common::ScalarType fileType = common::ScalarType::INTERNAL );

    /** 
     *  @brief Each processor writes an array of values at current file position with an individual offset
     */
    template<typename ValueType>
    void writeAll( const hmemo::HArray<ValueType>& local, const common::ScalarType fileType = common::ScalarType::INTERNAL );

    template<typename ValueType>
    void readSingle( ValueType val[], const IndexType n );

    template<typename ValueType>
    void readSingle( ValueType val[], const IndexType n, const common::ScalarType fileType );

    template<typename ValueType>
    void readSingle( ValueType& val, const common::ScalarType filetype = common::ScalarType::INTERNAL );

    template<typename ValueType>
    void readSingle( hmemo::HArray<ValueType>& val, const IndexType n, const common::ScalarType fileType = common::ScalarType::INTERNAL );

    template<typename ValueType>
    void readAll( ValueType val[], const IndexType n, const IndexType offset );

    template<typename ValueType>
    void readAll( ValueType local[], const IndexType size, const IndexType offset, const common::ScalarType fileType );

    template<typename ValueType>
    void readAll( hmemo::HArray<ValueType>& local, const IndexType size, const IndexType offset, 
                  const common::ScalarType stype = common::ScalarType::INTERNAL  );

    template<typename ValueType>
    void readAll( hmemo::HArray<ValueType>& local, const IndexType size, const common::ScalarType stype = common::ScalarType::INTERNAL );

    /**
     *  @brief Return the current pos in the file, is the number of bytes from beginning of the file.
     */
    inline size_t currentPos() const;

    inline const Communicator& getCommunicator() const;

    inline std::shared_ptr<const Communicator> getCommunicatorPtr() const;

protected:

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

    virtual void readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype ) = 0;

    virtual void readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype ) = 0;

    void set( const char* filename, size_t offset );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    CollectiveFile( CommunicatorPtr comm );

    std::shared_ptr<const Communicator> mComm;

    std::string mFileName;

    size_t mOffset;    // current file position
};

/* -------------------------------------------------------------------------- */
/*   Implementation of inline methods                                         */
/* -------------------------------------------------------------------------- */

size_t CollectiveFile::currentPos() const
{
    return mOffset;
}

std::shared_ptr<const Communicator> CollectiveFile::getCommunicatorPtr() const
{
    return mComm;
}

const Communicator& CollectiveFile::getCommunicator() const
{
    return *mComm;
}

}

}
