/**
 * @file FileIO.hpp
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
 * @brief Definition of abstract class for IO operations
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/Printable.hpp>

#include <scai/logging.hpp>

#include <string>

namespace scai
{

namespace common
{
    class Grid;
}

namespace lama
{

class _MatrixStorage;

/** This abstract base class specifies the methods that have to be
 *  implemented for read/write of matrix storage and vector arrays.
 *
 *  Each file format will derive from this class and provide the
 *  corresponding methods.
 *
 *  Note: For read and write operations the value type is required to
 *        make the correct output. But pure virtual methods can not
 *        have a template argument.
 *
 *        The struct IOWrapper is a helper struct that provides this
 *        funtionality.
 */

class COMMON_DLL_IMPORTEXPORT FileIO :

    public common::Factory<std::string, FileIO*>,
    public common::Printable

{

public:

    typedef enum
    {
        BINARY,             //!< forces the use of binary write
        FORMATTED,          //!< forces the use of formatted write
        DEFAULT_MODE        //!< keep it as it is set by default
    } FileMode;

    /** Constructor sets default values for mode */

    FileIO();

    virtual ~FileIO();

    /** Query if a certain file mode is supported
     *
     *  @param[in] mode either FORMATTED or BINARY
     *  @returns true if the file mode is supported
     */

    virtual bool isSupportedMode( const FileMode mode ) const = 0;

    /** Write matrix storage into a file.
     *
     *  @param[in] storage is the (local) matrix data that is written
     *  @param[in] fileName is the name of the output file
     *
     *  By default all data is written in exactly the same format as it
     *  is used in the storage. So for a storage with double values 8 bytes
     *  are needed for each entry, for a storage with float values 4 bytes.
     *
     *  - mBinary if true data is written binary
     *  - mScalarTypeIndex representation type used for row/column indexes (e.g. INT, LONG, ... ) within file
     *  - mScalarTypeData  representation type used for data within file
     */

    virtual void writeStorage( const _MatrixStorage& storage, const std::string& fileName ) = 0;

    /** Write 'dense' array of arbitrary type into a file.
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] fileName is the name of the output file
     *
     *  - mBinary if true data is written binary
     *  - mDataType output format used for non-zero values
     */
    virtual void writeArray( const hmemo::_HArray& array, const std::string& fileName ) = 0;

    /** Write multi-dimensional array of arbitrary type into a file with the shape
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] grid is the shape of the array
     *  @param[in] fileName is the name of the output file
     *
     *  - mBinary if true data is written binary
     *  - mDataType output format used for non-zero values
     */
    virtual void writeGridArray( const hmemo::_HArray& array, const common::Grid& grid, const std::string& fileName ) = 0;

    /** Write 'sparse' array of arbitrary type into a file.
     *
     *  @param[in] size is the full size of the array
     *  @param[in] indexes are the positions of the array with non-zero values
     *  @param[in] values are the values for the positions specified by indexes
     *  @param[in] fileName is the name of the output file
     *
     *  Note: the ZERO element of a sparse vector cannot be written into a file.
     */
    virtual void writeSparse( 
        const IndexType size, 
        const hmemo::HArray<IndexType>& indexes, 
        const hmemo::_HArray& values, 
        const std::string& fileName ) = 0;

    /** Get info about the storage stored in a file.
     *
     *  @param[out] numRows    number of rows for the storage in the file
     *  @param[out] numColumns number of columns for the storage in the file
     *  @param[out] numValues  number of non-zero values for the storage in the file
     *  @param[in]  fileName   number of non-zero values for the storage in the file
     */
    virtual void readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName ) = 0;

    /** Read (local) matrix storage from a file
     *
     *   - implicit conversion by reading from formatted file
     *   - for binary files the type must match the type of storage
     */

    /** Read in a matrix from a file but only a contiguous section of rows.
     *
     *  @param[out] storage  is the submatrix from the full matrix stored in the file
     *  @param[in]  fileName is the name of the input file containing the matrix.
     *  @param[in]  offsetRow index of first row to read
     *  @param[in]  nRows    number of rows to read
     *
     *  This routine can be used to read in one matrix from a single file with multiple processors
     *  where each processor reads in its local part.
     *
     *  This routine can also be used by a single processor that reads in the corresponding blocks and
     *  writes them to separate files.
     *
     *  The default implementation of the base class reads in the full storage and extracts the local
     *  part of it. Derived classes should implement solutions where it is not necessary to allocate memory
     *  for the full matrix but only for the corresponding block. In case of binary data, direct file access
     *  might be exploited to extract the needed data from the input file.
     */
    virtual void readStorage(
        _MatrixStorage& storage,
        const std::string& fileName,
        const IndexType offsetRow = 0,
        const IndexType nRows = invalidIndex ) = 0;

    /** Read in the size of an array saved in a file
     *
     *  @param[out] size     number of entries for the array saved in the file
     *  @param[in]  fileName C++ string containing the name of the file where the array is saved
     *  @throws common::Exception if file cannot be opened or if it does not contain an array
     */
    virtual void readArrayInfo( IndexType& size, const std::string& fileName ) = 0;

    /** Read in an array from a file in the corresponding format.
     *
     *  @param[out] array    container that will keep the array saved in file fileName
     *  @param[in]  fileName C++ string containing the name of the file where the array is saved
     *  @throws common::Exception if file cannot be opened or if it does not contain an array
     *
     *  If the value type of the array does not match the data stored in the file, an implicit
     *  type conversion is done.
     *
     *  If the file contains binary data, it is assumed that its type is the same as the value
     *  type of the array argument unless the environment variable ``SCAI_IO_TYPE`` has been set.
     */

    /** Read in a 'dense' array block from a file in the corresponding format.
     *
     *  @param[out] array    will contain the corresponding array values
     *  @param[in]  fileName name of the input file with array data
     *  @param[in]  offset   first entry to read
     *  @param[in]  n        number of entries to read, invalidIndex stands for all remaining entries
     *
     *  This method has exactly the same behavior as readArray but with the difference that only
     *  a part of the array is read.
     */
    virtual void readArray(
        hmemo::_HArray& array,
        const std::string& fileName,
        const IndexType offset = 0,
        const IndexType n = invalidIndex ) = 0;

    virtual void readGridArray(
        hmemo::_HArray& array,
        common::Grid& grid,
        const std::string& fileName ) = 0;

    /** Read in a 'sparse' array from a file in the corresponding format.
     *
     *  @param[out] size is the size of the array
     *  @param[out] indexes are the positions with non-zero values
     *  @param[out] values are the values at the corresponding positions.
     *  @param[in]  fileName of the input file with array data.
     *
     *  This method must be implemented by each derived class. 
     */
    virtual void readSparse( 
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values,
        const std::string& fileName ) = 0;

    /** File suffix can be used to decide about choice of FileIO class */

    virtual std::string getMatrixFileSuffix() const = 0;

    virtual std::string getVectorFileSuffix() const = 0;

    /** Help routine to delete file and maybe joint files
     *
     *  @param[in] fileName file to delete
     *  @returns   0 on success
     */
    virtual int deleteFile( const std::string& fileName );

    /** Setter for representation type used for indexes in file. */

    void setIndexType( common::ScalarType type );

    /** Setter for representation type used for data values in file. */

    void setDataType( common::ScalarType type );

    /** Enable/Disable binary mode. */

    void setMode( FileMode mode );

    /** Enable/Disable append mode.
     *
     *  For the append mode, write operations will append data in existing files
     */
    void enableAppendMode( bool flag );

    /** Help routine that determines the availability of a given file by its name. */

    static bool fileExists( const std::string& fileName );

    /** Query if a file has a certain suffix */

    static bool hasSuffix( const std::string& fileName, const std::string& suffix );

    /** Return the suffix of a filename, with point, e.g. ".mtx" for "a.mtx" */

    static std::string getSuffix( const std::string& fileName );

    /** Help routine to remove file and maybe joint files
     *
     *  @param[in] fileName file to delete
     *  @returns   0 on success
     */
    static int removeFile( const std::string& fileName );

    /** Static method to write a dense array into a file.
     *
     *  @param[in] array     is the array that is saved
     *  @param[in] fileName  is the name of output file, suffix decides about Handler
     *  @param[in] dataType  specifies the type to be used for representation in output file
     *
     *  If the optional argument dataType is not set, array.getValueType() is used.
     */
    static void write(
        const hmemo::_HArray& array,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    /** Static method to write a sparse array into a file.
     *
     *  @param[in] size, values, indexes represent the sparse array
     *  @param[in] fileName  is the name of output file, suffix decides about Handler
     *  @param[in] dataType  specifies the type to be used for representation in output file
     *
     *  If the optional argument dataType is not set, array.getValueType() is used.
     */
    static void write(
        const IndexType size,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::_HArray& values,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    /** Static method to read an array or a contiguous section of an array from a file.
     *
     *  @param[out] array will contain the array data
     *  @param[in]  fileName is the name of the input file where the array is saved
     *  @param[in]  dataType specifies the type that has been used for the values in the input file
     *  @param[in]  first index of the first element to read, defaults to 0
     *  @param[in]  n     number of elements, default is invalidIndex that stands for up to the end
     *
     *  \code
     *      HArray<double> data;
     *      FileIO::read( data, "myData.txt" )           // reads the full array
     *      FileIO::read( data, "myData.txt", 50 )       // reads array but skips first 50 entries
     *      FileIO::read( data, "myData.txt", 50, 10 )   // reads for pos 50 next 10 elements
     *  \endcode
     */
    static void read(
        hmemo::_HArray& array,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL,
        const IndexType first = 0,
        const IndexType n = invalidIndex );

    static void read(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& array,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    /** Stati method to read a multi-dimensional array from a file. */

    static void read(
        hmemo::_HArray& array,
        common::Grid& grid,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    /** Return the size of an array saved in a file. */

    static IndexType getArraySize( const std::string& fileName );

    /** Return the size of a storage saved in a file. */

    static IndexType getStorageSize( const std::string& fileName );

    /** Override Printable::writeAt */

    void writeAt( std::ostream& stream ) const;

protected:

    /** write the global settings in the stream, useful for derived classes */

    void writeMode( std::ostream& stream ) const;

    /** get precision used for formatted output of values
     *
     *  @param[in] valueType is the actual type of values to be printed
     */

    int getDataPrecision( common::ScalarType valueType );

    FileMode mFileMode;                     //!< if true binary mode must be respected
    bool mAppendMode;                       //!< if true output is appended to existing files

    common::ScalarType mScalarTypeIndex; //!< representation type of row indexes
    common::ScalarType mScalarTypeData;  //!< representation type of array/matrix value entries

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
