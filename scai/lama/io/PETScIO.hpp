/**
 * @file PETScIO.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Structure that contains IO routines for PETSC
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama/io/FileIO.hpp>
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

    /** Implementation of pure virtual method FileIO::writeStorage */

    void writeStorage( const _MatrixStorage& storage, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage(
        _MatrixStorage& storage,
        const std::string& fileName,
        const IndexType offsetRow,
        const IndexType nRows );

    /** Implementation of FileIO::writeArray */

    void writeArray( const hmemo::_HArray& array, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::writeSparse  */

    virtual void writeSparse(
        const IndexType size,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::_HArray& array,
        const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readArray using same defaults */

    virtual void readArray(
        hmemo::_HArray& array,
        const std::string& fileName,
        const IndexType offset = 0,
        const IndexType n = invalidIndex );

    /** Implementation of pure virtual method FileIO::readSparse */

    virtual void readSparse(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values,
        const std::string& fileName );

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

    /** Implementation of pure methdod FileIO::readStorageInfo */

    virtual void readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName );

    /** Implementation of pure methdod FileIO::readArrayInfo */

    virtual void readArrayInfo( IndexType& size, const std::string& fileName );

    /** Implementation of writing array with grid information */

    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName );

    void readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& outputFileName );

public:

    /** Typed version of writeStorage called via IOWrapper */

    template<typename ValueType>
    void writeStorageImpl( const MatrixStorage<ValueType>& storage, const std::string& fileName );

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName, const IndexType firstRow, const IndexType nRows );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName );

    /** Typed version of the writeSparse */

    template<typename ValueType>
    void writeSparseImpl(
        const IndexType size,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::HArray<ValueType>& values,
        const std::string& fileName );

    /** Typed version of readArray */

    template<typename ValueType>
    void readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName, const IndexType first, const IndexType n );

    /** Typed version of readSparse */

    template<typename ValueType>
    void readSparseImpl(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<ValueType>& values,
        const std::string& fileName );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

}

}
