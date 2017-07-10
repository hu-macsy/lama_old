/**
 * @file BitmapIO.hpp
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

    /** Implementation of pure method FileIO::readStorageInfo */

    virtual void readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName );

    /** Implementation of pure methdod FileIO::readArrayInfo */

    virtual void readArrayInfo( IndexType& size, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::writeStorage for all derived classes */

    virtual void writeStorage( const _MatrixStorage& storage, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readStorage  */

    virtual void readStorage( _MatrixStorage& storage, const std::string& fileName, const IndexType offsetRow, const IndexType nRows );

    /** Implementation of pure virtual method FileIO::writeArray  */

    virtual void writeArray( const hmemo::_HArray& array, const std::string& fileName );

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
        const IndexType n = nIndex );

    /** Implementation of pure virtual method FileIO::readSparse */

    virtual void readSparse(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values,
        const std::string& fileName );

    /** Default implementation for query matrix file suffix, is createValue of derived class */

    virtual std::string getMatrixFileSuffix() const;

    /** Default implementation for query vector file suffix, is createValue of derived class */

    virtual std::string getVectorFileSuffix() const;

    /** Implementation of pure methdod FileIO::isSupportedMode */

    virtual bool isSupportedMode( const FileMode mode ) const;

    void readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& outputFileName );

    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName );

    /** Typed version of BitmapIO::read */

    template<typename ValueType>
    void readGridImpl( hmemo::HArray<ValueType>& data, common::Grid& grid, const std::string& outputFileName );

    /** Typed version of BitmapIO::write */

    template<typename ValueType>
    void writeGridImpl( const hmemo::HArray<ValueType>& data, const common::Grid& grid, const std::string& outputFileName );

    /** Implementation for Printable.:writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    // registration key for factory

    static std::string createValue();

    // static method to create an FileIO object for this derived class

    static FileIO* create();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

} /* end namespace lama */

} /* end namespace scai */
