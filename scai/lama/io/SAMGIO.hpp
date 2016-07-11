/**
 * @file SAMGIO.hpp
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
 * @brief Dervied FileIO class that implements IO routines for the SAMG file format. 
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama/io/CRTPFileIO.hpp>

namespace scai
{

namespace lama
{

/** The SAMG format supports binary and formatted IO.
 *
 *  The SAMG format uses two files, one for header information and one with
 *  the data itself. Only the header file name has to be specified for
 *  read/write operations.
 */

class SAMGIO : 

    public CRTPFileIO<SAMGIO>,         // use type conversions
    public FileIO::Register<SAMGIO>    // register at factory
{

public:

    /** Constructor might reset default values */

    SAMGIO();

    /** Implementation of pure methdod FileIO::isSupportedMode */

    virtual bool isSupportedMode( const FileMode mode ) const;

    /** File suffix is used to decide about choice of output class */

    virtual std::string getVectorFileSuffix() const;

    virtual std::string getMatrixFileSuffix() const;

    /** Implementation for Printable.:writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    // static method to create an FileIO object for this derived class

    static FileIO* create();

    // registration key for factory

    static std::string createValue();

    /** Override default implementation as also data files should be deleted. */

    virtual int deleteFile( const std::string& fileName );

public:
 
    /** Typed version of writeStorage
     *
     *  This method must be available for implementation of
     *  CRTPFileIO::writeStorage
     */

    template<typename ValueType>
    void writeStorageImpl( const MatrixStorage<ValueType>& storage, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of readArray */

    template<typename ValueType>
    void readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName )
    __attribute( ( noinline ) );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class

private:

    /** Guard class for an additional registration with the vector file suffix. */

    class Guard
    {
    public:

        Guard();
        ~Guard();
    };

    static Guard mGuard;

    /** Own routines for read/write of the header file, no template parameters required */

    void readMatrixHeader( IndexType& numRows, IndexType& numValues, bool& binary, const std::string& fileName );

    void writeMatrixHeader( const IndexType numRows, const IndexType numValues, const bool binary, const std::string& fileName );

    void readVectorHeader( IndexType& n, IndexType& typeSize, bool& binary, const std::string& fileName );

    void writeVectorHeader( const IndexType n, const IndexType typeSize, const bool binary, const std::string& fileName );

};

}

}
