/**
 * @file TextIO.hpp
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
 * @brief Structure that contains IO routines for simple ASCII text files
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

/** This file format stores just the COO data of a matrix.
 *
 *   - there is not header at all
 *   - number of non-zero matrix entries is given by number of lines
 *   - size of matrix is given by maximal values for row and for column indexes
 *   - size of vector is just given by number of lines
 *
 *   It is very useful to read dumped matrices of Matlab.
 *
 *   /code
 *   data_mat = [ia ja,real(val),imag(val)];
 *   save -ascii dataMatrix.txt data_mat;
 *
 *   data_rhs = [real(b),imag(b)];
 *   save -ascii datRHS.txt data_rhs;
 *   /endcode
 *
 *   /code
 *   solver.exe dataMatrix.txt datRHS.txt
 *   /endcode
 */

class TextIO :

    public FileIO,    
    public FileIO::Register<TextIO>    // register at factory
{

public:

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

public:

    /** Typed version of writeStorage */

    template<typename ValueType>
    void writeStorageImpl( const MatrixStorage<ValueType>& storage, const std::string& fileName );

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName, const IndexType firstRow, const IndexType nRows );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName );

    /** Typed version of writeArray for sparse arrays. */

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

    /** Implementation of writing array with grid information */

    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName );

    void readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& outputFileName );

private:

    /** Help routine to read storage data */

    template<typename ValueType>
    void readData(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>* vals,
        const IndexType nnz,
        const std::string& fileName );

    /** Method to count number of lines of a text file and the maximal number of entries in one line
     *
     *  @param[out]  nLines will will contain the number of lines the file has
     *  @param[out]  nEntries is maximal number of entries in one line
     *  @param[in]   fileName is the name of the file
     *
     *  Note: it might be possible that one line contains less than 'nEntries' entries
     */
    void checkTextFile( IndexType& nLines, IndexType& nEntries, const char* fileName );
};

}

}
