/**
 * @file FileIO.hpp
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
 *        The class CRTPFileIO is a helper class that provides this 
 *        funtionality.
 */

class FileIO : 

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

    /** Write array of arbitrary type into a file.
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] fileName is the name of the output file
     *
     *  - mBinary if true data is written binary
     *  - mDataType output format used for non-zero values
     */

    virtual void writeArray( const hmemo::_HArray& array, const std::string& fileName ) = 0;

    /** Get info about the storage stored in a file. 
     *
     *  @param[out] numRows    number of rows for the storage in the file
     *  @param[out] numColumns number of columns for the storage in the file
     *  @param[out] numValues  number of non-zero values for the storage in the file
     */
    virtual void readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName ) = 0;

    /** Read (local) matrix storage from a file
     *
     *   - implicit conversion by reading from formatted file
     *   - for binary files the type must match the type of storage
     */

    virtual void readStorage( _MatrixStorage& storage, const std::string& fileName ) = 0;

    /** Read in a matrix from a file but only a contiguous section of rows. 
     *
     *  @param[out] storage  is the submatrix from the full matrix stored in the file
     *  @param[in]  fileName is the name of the input file containing the matrix.
     *  @param[in]  firstRow index of first row to read
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
    virtual void readStorageBlock( 
        _MatrixStorage& storage, 
        const std::string& fileName,
        const IndexType firstRow,
        const IndexType nRows ) = 0;

    virtual void readArray( hmemo::_HArray& array, const std::string& fileName ) = 0; 

    virtual void readArrayBlock( 
        hmemo::_HArray& array,
        const std::string& fileName,
        const IndexType first,
        const IndexType n ) = 0;

    /** File suffix can be used to decide about choice of FileIO class */

    virtual std::string getMatrixFileSuffix() const = 0;

    virtual std::string getVectorFileSuffix() const = 0;

    /** Help routine to delete file and maybe joint files 
     *
     *  @param[in] fileName file to delete
     *  @returns   0 on success
     */
    virtual int deleteFile( const std::string& fileName ) = 0;

    /** Setter for representation type used for indexes in file. */

    void setIndexType( common::scalar::ScalarType type );

    /** Setter for representation type used for data values in file. */

    void setDataType( common::scalar::ScalarType type );

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

    static void write( const hmemo::_HArray& array, const std::string& fileName );

    static void read( hmemo::_HArray& array, const std::string& fileName );

    /** Override Printable::writeAt */

    void writeAt( std::ostream& stream ) const;

protected:

    /** write the global settings in the stream, useful for derived classes */

    void writeMode( std::ostream& stream ) const;

    /** get precision used for formatted output of values
     *
     *  @param[in] valueType is the actual type of values to be printed
     */

    int getDataPrecision( common::scalar::ScalarType valueType );

    FileMode mFileMode;                     //!< if true binary mode must be respected
    bool mAppendMode;                       //!< if true output is appended to existing files

    common::scalar::ScalarType mScalarTypeIndex; //!< representation type of row indexes
    common::scalar::ScalarType mScalarTypeData;  //!< representation type of array/matrix value entries

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
