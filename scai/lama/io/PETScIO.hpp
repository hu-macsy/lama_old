/**
 * @file PETScIO.hpp
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
 * @brief Structure that contains IO routines for PETSC
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

namespace scai
{

namespace lama
{

/** This file format supports the binary format used by PetSC.
 *
 *   - header information is just at the beginning of the file
 *   - uses CSR format, but the sizes array and not the offsets
 *   - stores data always in BIG endian (x86 has LITTLE endian)
 */

class PETScIO :

    public FileIO, 
    public FileIO::Register<PETScIO>    // register at factory

{

public:

    /** Constructor might reset default values */

    PETScIO();

    /** Implementation of pure virtual method FileIO::open */

    virtual void open( const char* fileName, const char* fileMode );

    /** Implementation of pure virtual method FileIO::close */

    virtual void close();

    /** Implementation of pure virtual method FileIO::writeStorage */

    void writeStorage( const _MatrixStorage& storage );

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage( _MatrixStorage& storage );

    /** Implementation of FileIO::writeArray */

    void writeArray( const hmemo::_HArray& array );

    /** Implementation of pure virtual method FileIO::writeSparse  */

    virtual void writeSparse(
        const IndexType size,
        const void* zero,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::_HArray& array );

    /** Implementation of pure virtual method FileIO::readArray using same defaults */

    virtual void readArray( hmemo::_HArray& array );

    /** Implementation of pure virtual method FileIO::readSparse */

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

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage );

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

private:

    class IOStream mFile;    // used file

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

}

}
