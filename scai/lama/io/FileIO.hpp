/**
 * @file FileIO.hpp
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

/**
 *  @brief Enumeration class for different file modes
 */
enum class FileMode
{
    BINARY,             //!< binary
    FORMATTED,          //!< formatted text
    DEFAULT             //!< keep it as it is set by default
};

/**
 *  Enum type for entities that might be stored in a file.
 */
enum class IOItem
{
    DENSE_DATA,         //!< dense vector
    GRID_DATA,          //!< grid vector
    SPARSE_DATA,        //!< sparse vector
    CSR_MATRIX,         //!< sparse matrix in CSR format
    COO_MATRIX,         //!< sparse matrix in COO format
    DENSE_MATRIX,       //!< sparse matrix in dense format
    UNKNOWN             //!< unknown entry, i.e. any other stuff
};

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

    /** Constructor sets default values for mode */

    FileIO();

    virtual ~FileIO();

    /**
     *  @brief Pure method to open a file.
     *
     *  @param[in] fileName specifies the name of the file
     *  @param[in] fileMode is either "r" for reading a file, "w" for writing into it, or "a" for append to it
     *
     *  This method might throw an IOException if the file could not be opened or if a mode is not supported.
     */
    virtual void open( const char* fileName, const char* fileMode ) = 0;

    /**
     *  @brief Pure method to close the currently open file.
     *
     *  This method might throw an IOEception if no file is open at all.
     */
    virtual void close() = 0;

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

    virtual void writeStorage( const _MatrixStorage& storage ) = 0;

    /** Write 'dense' array of arbitrary type into a file.
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] fileName is the name of the output file
     *
     *  - mBinary if true data is written binary
     *  - mDataType output format used for non-zero values
     */
    virtual void writeArray( const hmemo::_HArray& array ) = 0;

    /** Write multi-dimensional array of arbitrary type into a file with the shape
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] grid is the shape of the array
     *  @param[in] fileName is the name of the output file
     *
     *  - mBinary if true data is written binary
     *  - mDataType output format used for non-zero values
     */
    virtual void writeGridArray( const hmemo::_HArray& array, const common::Grid& grid ) = 0;

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
        const hmemo::_HArray& values ) = 0;

    /** Get info about the storage stored in a file.
     *
     *  @param[out] numRows    number of rows for the storage in the file
     *  @param[out] numColumns number of columns for the storage in the file
     *  @param[out] numValues  number of non-zero values for the storage in the file
     *
     */
    virtual void getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues ) = 0;

    /** Read in a matrix storage from the opened file
     *
     *  @param[out] storage  is the submatrix from the full matrix stored in the file
     */
    virtual void readStorage( _MatrixStorage& storage ) = 0;

    /** Read in the size of an array saved in a file
     *
     *  @param[out] size     number of entries for the array saved in the file
     *  @throws common::Exception if file cannot be opened or if it does not contain an array
     */
    virtual void getArrayInfo( IndexType& size ) = 0;

    /** Read in a 'dense' array block from a file in the corresponding format.
     *
     *  @param[out] array    will contain the corresponding array values
     */
    virtual void readArray( hmemo::_HArray& array ) = 0;

    virtual void readGridArray(
        hmemo::_HArray& array,
        common::Grid& grid ) = 0;

    /** Read in a 'sparse' array from a file in the corresponding format.
     *
     *  @param[out] size is the size of the array
     *  @param[out] indexes are the positions with non-zero values
     *  @param[out] values are the values at the corresponding positions.
     *
     *  This method must be implemented by each derived class. 
     */
    virtual void readSparse( 
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values ) = 0;

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

    /** Static method to read an array from a file
     *
     *  @param[out] array will contain the array data
     *  @param[in]  fileName is the name of the input file where the array is saved
     *  @param[in]  dataType specifies the type that has been used for the values in the input file
     *
     *  \code
     *      HArray<double> data;
     *      FileIO::read( data, "myData.txt" )           // reads the full array
     *  \endcode
     */
    static void read(
        hmemo::_HArray& array,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    static void read(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& array,
        const std::string& fileName,
        const common::ScalarType dataType = common::ScalarType::INTERNAL );

    /** Static method to read a multi-dimensional array from a file. */

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
