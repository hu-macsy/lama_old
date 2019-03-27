/**
 * @file lama/io/BitmapIO.hpp
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
 * @brief Definition of routines to read/write image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/logging.hpp>
#include <scai/lama/io/FileIO.hpp>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT BitmapIO :

    public FileIO,
    public FileIO::Register<BitmapIO>    // register at factory

{

public:

    /** Implementation of pure virtual method FileIO::openIt */

    virtual void openIt( const std::string& fileName, const char* fileMode );

    /** Implementation of pure virtual method FileIO::closeIt */

    virtual void closeIt();

    /** Implementation of pure method FileIO::getStorageInfo */

    virtual void getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues );

    /** Implementation of pure methdod FileIO::getArrayInfo */

    virtual void getArrayInfo( IndexType& size );

    /** Implementation of pure virtual method FileIO::writeStorage for all derived classes */

    virtual void writeStorage( const _MatrixStorage& storage );

    /** Implementation of pure virtual method FileIO::readStorage  */

    virtual void readStorage( _MatrixStorage& storage );

    /** Implementation of pure virtual method FileIO::writeArray  */

    virtual void writeArray( const hmemo::_HArray& array );

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

    /** Default implementation for query matrix file suffix, is createValue of derived class */

    virtual std::string getMatrixFileSuffix() const;

    /** Default implementation for query vector file suffix, is createValue of derived class */

    virtual std::string getVectorFileSuffix() const;

    /** Implementation of pure methdod FileIO::isSupportedMode */

    virtual bool isSupportedMode( const FileMode mode ) const;

    /** Override methdod FileIO::readGridArray */

    void readGridArray( hmemo::_HArray& data, common::Grid& grid );

    /** Override methdod FileIO::writeGridArray */

    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid );

    /** Typed version of BitmapIO::read */

    template<typename ValueType>
    void readGridImpl( hmemo::HArray<ValueType>& data, common::Grid& grid );

    /** Typed version of BitmapIO::write */

    template<typename ValueType>
    void writeGridImpl( const hmemo::HArray<ValueType>& data, const common::Grid& grid );

    /** Implementation for Printable.:writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    // registration key for factory

    static std::string createValue();

    // static method to create an FileIO object for this derived class

    static FileIO* create();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    FILE* mFile;

};

} /* end namespace lama */

} /* end namespace scai */
