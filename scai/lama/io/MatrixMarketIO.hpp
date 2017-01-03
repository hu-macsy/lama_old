/**
 * @file MatrixMarketIO.hpp
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
 * @brief Structure that contains IO routines for PETSC
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama/io/CRTPFileIO.hpp>

namespace scai
{

namespace lama
{

class MatrixMarketIO :

    public CRTPFileIO<MatrixMarketIO>,         // use type conversions
    public FileIO::Register<MatrixMarketIO>    // register at factory
{

public:

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

    /** Implementation of pure methdod FileIO::readStorageInfo */

    virtual void readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName );

    /** Implementation of pure methdod FileIO::readArrayInfo */

    virtual void readArrayInfo( IndexType& size, const std::string& fileName );

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
    void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName, const IndexType firstRow, const IndexType nRows )
    __attribute( ( noinline ) );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of readArray */

    template<typename ValueType>
    void readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName, const IndexType first, const IndexType n )
    __attribute( ( noinline ) );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
 
private:

    /** Enumeration type for the different symmetry flags in the Matrix Market file */

    typedef enum
    {
        GENERAL,
        SYMMETRIC,
        HERMITIAN,
        SKEW_SYMMETRIC
    } Symmetry;

    /** Conversion of enum value to string */

    const char* symmetry2str( const Symmetry symmetry );

    void writeMMHeader(
        class IOStream& outFile,
        const bool vector,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const Symmetry symmetry,
        const common::scalar::ScalarType dataType );

    void readMMHeader(
        IndexType& numRows,
        IndexType& numColumns,
        IndexType& numValues,
        common::scalar::ScalarType& dataType,
        Symmetry& symmetry,
        class IOStream& inFile );

    /** Extend COO data by adding all symmetric data */

    template<typename ValueType>
    void addSymmetricEntries(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& vals, 
        bool conjFlag );

    /** Check for symmetry */
    template<typename ValueType>
    Symmetry checkSymmetry( 
        const hmemo::HArray<IndexType>& cooIA, 
        const hmemo::HArray<IndexType>& cooJA, 
        const hmemo::HArray<ValueType>& cooValues );

    template<typename ValueType>
    void readMMArray(
        class IOStream& inFile,
        hmemo::HArray<ValueType>& data,
        const IndexType numRows,
        const IndexType numColumns,
        const Symmetry symmetry );
};

}

}
