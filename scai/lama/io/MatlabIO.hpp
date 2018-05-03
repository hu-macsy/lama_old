/**
 * @file MatlabIO.hpp
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
 * @brief Structure that contains IO routines for MAT-File Format
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

/** This file I/O class supports the Level 5 MAT-File Format of MATLAB.
 *
 *  It allows to exchange data between MATLAB and LAMA applications.
 *
 *   /code
 *   save example.mat matrix
 *   load example.mat
 *   /endcode
 *
 *   /code
 *   CSRStorge<ValueType> storage;
 *   storage.readFromFile( "example.mat" );
 *   storage.writeToFile( "example.mat" );
 *   /endcode
 *
 *  The following topics should be observed:
 *
 *  - LAMA can only read and write one single data element from a MATLAB file.
 *  - The name of the data element is is ignored when reading the element and
 *    each written element gets the name "LAMA".
 *  - As the data type is stored for each element in the file, the SCAI_IO_TYPE
 *    is ignored, i.e. each array/storage is written exactly in the format it is
 *    and there might be an implicit type conversion during the read.
 *  - Only the SCAI_IO_TYPE=PATTERN will be handled and in this case no sparse matrix
 *    values are written, only the row and column indexes
 *  - The data types LongDouble and ComplexLongDouble are not really supported
 *    by the Level 5 MAT-File format but are full supported here by using a reserved
 *    value of the MAT-File Data Types.
 */

class MatlabIO :

    public FileIO,        
    public FileIO::Register<MatlabIO>    // register at factory
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

    /** Implementation of pure method FileIO::writeGridArray.
     *
     *  The MATLAB format supports writing the shape information.
     */  
    void writeGridArray( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName );

    /** Typed version of MatlabIO::writeGridArray */

    template<typename ValueType>
    void writeGridImpl( const hmemo::HArray<ValueType>& data, const common::Grid& grid, const std::string& outputFileName );

    void readGridArray( hmemo::_HArray& data, common::Grid& grid, const std::string& outputFileName );

    /** Typed version of MatlabIO::readGridArray */

    template<typename ValueType>
    void readGridImpl( hmemo::HArray<ValueType>& data, common::Grid& grid, const std::string& outputFileName );

private:

    template <typename ValueType>
    uint32_t getSparseStorage( MatrixStorage<ValueType>& storage,
                               const IndexType dims[2], const IndexType nnz,
                               bool isComplex,
                               const char* dataElementPtr, uint32_t nBytes );


    template <typename ValueType>
    uint32_t getStructStorage( MatrixStorage<ValueType>& storage, const char* dataElementPtr, uint32_t nBytes );

    /** Help routine to read storage data, uses IO types */

    template <typename ValueType>
    void getStorage( MatrixStorage<ValueType>& storage, const char* dataElementPtr, uint32_t nBytes );

    /** Help routine to read a heterogeneous array from input data
     *
     *  @param[in] data is the input data with header info about size and type
     *  @param[in] len  is the size of the input data, used for checks to avoid out-of-range
     *  @param[out] array will contain the read data, might be converted to ValueType
     */
    template<typename ValueType>
    static uint32_t getArrayData( hmemo::HArray<ValueType>& array, const char* data, uint32_t len );

    template<typename ValueType>
    static void readMATArray( hmemo::HArray<ValueType>& array, const char* data, const uint32_t mxType, const uint32_t nbytes );

    /** Help routine to read a heterogeneous array
     *
     *  @tparam ArrayType is the value type of the heterogeneous array
     *  @tparam DataType is the type of data from the file
     */
    template<typename ArrayType, typename DataType>
    static void readMATArrayImpl( hmemo::HArray<ArrayType>& array, const void* data, IndexType nBytes );

    /** Write a heterogeneous array as MAT format into output file. */

    template<typename ValueType>
    static uint32_t writeArrayData( class MATIOStream& outFile, const hmemo::HArray<ValueType>& array, bool dryRun );

    /** Write dense array or matrix as complete element in a Matlab file. */

    template<typename ValueType>
    void writeDenseArray( class MATIOStream& outFile, const hmemo::HArray<ValueType>& array, IndexType dims[] );

    template<typename ValueType>
    void writeDenseGrid( MATIOStream& outFile, const hmemo::HArray<ValueType>& array, const common::Grid& grid );

};

}

}
