/**
 * @file MatrixMarketIO.hpp
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

namespace scai
{

namespace lama
{

template<typename ValueType>
class DenseStorage;

template<typename ValueType>
class MatrixStorage;

/** Enumeration type for the different symmetry flags in the _Matrix Market file */

enum class Symmetry
{
    GENERAL,               //!< no symmetry used
    SYMMETRIC,             //!< a( i, j ) and a( j, i ) are always same
    HERMITIAN,             //!< a( i, j ) == conj( a( i, j  ) )
    SKEW_SYMMETRIC         //!< not exploited here
};

/** 
 *  @brief Structure that contains all data available in the header of a matrix market file
 */
struct MMHeader
{
    common::ScalarType mmType;   //!< specifies the type of the data, e.g. real or complex

    bool      isVector;          //!< either vector (1D) or matrix (2D)

    IndexType mNumRows;           //!< number of rows
    IndexType mNumColumns;        //!< number of columns, 1 if isVector
    IndexType mNumValues;         //!< number of values if coordinates are used (sparse format)

    Symmetry symmetry;           //!< only used for matrix, should be GENERAL for vector

    /** Constructor of a header for a dense vector */

    MMHeader( const common::ScalarType dataType, const IndexType size );

    /** Constructor of a header for a dense matrix */

    MMHeader( const common::ScalarType dataType, const IndexType numRows, const IndexType numColumns );
};

class MatrixMarketIO :

    public FileIO,
    public FileIO::Register<MatrixMarketIO>    // register at factory
{

public:

    /** Implementation of pure virtual method FileIO::openIt */

    virtual void openIt( const std::string& fileName, const char* fileMode );

    /** Implementation of pure virtual method FileIO::close */

    virtual void closeIt();

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

    /** Implementation of FileIO::readGridArray */

    void readGridArray( hmemo::_HArray& data, common::Grid& grid );

public:

    /** Typed version of writeStorage
     *
     *  This method must be available for implementation of
     */

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

    /** Typed version of readGridArray */

    template<typename ValueType>
    void readGridImpl( hmemo::HArray<ValueType>& data, common::Grid& grid );

    /** Typed version of writeGridArray */

    template<typename ValueType>
    void writeGridImpl( const hmemo::HArray<ValueType>& data, const common::Grid& grid );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class

private:

    class IOStream mFile;    // used file

    /** Conversion of enum value to string */

    const char* symmetry2str( const Symmetry symmetry );

    void writeMMHeader( class IOStream& outFile, const MMHeader& header );

    /** 
     *  @brief read the complete header of a matrix-market file
     */
    MMHeader readMMHeader();

    /** 
     *  @brief same as readMMHeader but keep the current file position
     */
    MMHeader getMMHeader();

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
        hmemo::HArray<ValueType>& data,
        const IndexType numRows,
        const IndexType numColumns,
        const Symmetry symmetry );

    template<typename ValueType>
    void readVectorCoordinates(
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<ValueType>& values,
        const IndexType numValues,
        const bool isVector,
        common::ScalarType mmType );

    template<typename ValueType>
    void writeDenseMatrix( const DenseStorage<ValueType>& storage );

    template<typename ValueType>
    void writeArray2D( const IndexType numRows, const IndexType numColumns, const hmemo::HArray<ValueType>& data );

};

}

}
