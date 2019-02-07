/**
 * @file LamaIO.hpp
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
 * @brief Proprietray binary LAMA IO format
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/IOStream.hpp>

#include <scai/dmemo/CollectiveFile.hpp>

namespace scai
{

namespace lama
{

// forward declarations

class _MatrixStorage;
template<typename ValueType> class MatrixStorage;
template<typename ValueType> class CSRStorage;

/** 
 *  Collective/Parallel I/O format used by LAMA
 *
 *   - supports collective mode, i.e. parallel read/write at different offsets
 *   - supports dense vector  
 *   - supports sparse vector
 *   - supports dense matrix
 *   - supports csr matrix
 *   - supports coo matrix
 *
 *  Note: coo/sparse format allows to write vector/matrices that are not 
 *        block distributed (global coordinates)
 *        dense vector/matrix and csr must be block distributed
 */
class LamaIO :

    public FileIO, 
    public FileIO::Register<LamaIO>    // register at factory

{

public:

    /** Constructor might reset default values */

    LamaIO();

    /** Implementation of pure virtual method FileIO::openIt */

    virtual void openIt( const std::string& fileName, const char* fileMode );

    /** Implementation of pure virtual method FileIO::close */

    virtual void closeIt();

    /** Override default method as this class fully supports collective I/O */

    virtual bool hasCollectiveIO() const;

    /** Implementation of pure virtual method FileIO::writeStorage */

    void writeStorage( const _MatrixStorage& storage );

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage( _MatrixStorage& storage );

    /** Implementation of FileIO::writeArray */

    void writeArray( const hmemo::_HArray& array );

    /** Override virtual method FileIO::writeSparse as it is supported */

    virtual void writeSparse(
        const IndexType size,
        const void* zero,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::_HArray& array );

    /** Implementation of pure virtual method FileIO::readArray using same defaults */

    virtual void readArray( hmemo::_HArray& array );

    /** Override virtual method FileIO::readSparse as it is supported */

    virtual void readSparse(
        IndexType& size,
        void* zero,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values );

    /** Implementation of FileIO::getMatrixFileSuffix */

    std::string getMatrixFileSuffix() const;

    /** Implementation of FileIO::getVectorFileSuffix */

    std::string getVectorFileSuffix() const;

    /** Implementation of pure methdod FileIO::isSupportedMode */

    virtual bool isSupportedMode( const FileMode mode ) const;

    /** Implementation for Printable.:writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    // static method to create an FileIO object for this derived class

    static FileIO* create();

    // registration key for factory

    static std::string createValue();

    /** Implementation of pure methdod FileIO::getStorageInfo */

    virtual void getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues );

    /** Implementation of pure methdod FileIO::getArrayInfo */

    virtual void getArrayInfo( IndexType& size );

    /** Implementation of writing array with grid information */

    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid );

    void readGridArray( hmemo::_HArray& data, common::Grid& grid );

    /** Typed version of writeStorage called via IOWrapper */

    template<typename ValueType>
    void writeStorageImpl( const MatrixStorage<ValueType>& storage );

    template<typename ValueType>
    void writeCSR( const CSRStorage<ValueType>& storage );

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage );

    template<typename ValueType>
    void readCSR( CSRStorage<ValueType>& storage );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array );

    /** Typed version of the writeSparse */

    template<typename ValueType>
    void writeSparseImpl(
        const IndexType size,
        const ValueType& zero,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::HArray<ValueType>& values );

    /** Typed version of readArray */

    template<typename ValueType>
    void readArrayImpl( hmemo::HArray<ValueType>& array );

    /** Typed version of readSparse */

    template<typename ValueType>
    void readSparseImpl(
        IndexType& size,
        ValueType& zero,
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<ValueType>& values );

    void writeHeader( const int classId );

private:

    std::unique_ptr<class dmemo::CollectiveFile> mFile;    // used file

    std::string mFileName;

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

}

}
