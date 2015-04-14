/**
 * @file StorageIO.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Class providing IO routines for matrix storage
 * @author Thomas Brandes
 * @date 31.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_STORAGE_IO_HPP_
#define LAMA_STORAGE_IO_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMAArray.hpp>
#include <lama/io/FileType.hpp>

#include <fstream>

namespace lama
{

/* -------------------------------------------------------------------------- */

/** This class provides static utility methods for splitting matrix storage into a
 *  local and a halo part.
 *
 *  Due to a column distribution the storage is divided into a local part (having
 *  the local columns) and a halo part (for the non-local columns). Furthermore,
 *  it builds the halo for exchanging the non-local values between processors.
 */
class LAMA_DLL_IMPORTEXPORT _StorageIO
{
public:

    /** This method writes a header file for a CSR storage.
     *
     *  @param[in] numRows    is number of rows of the storage
     *  @param[in] numValues  is number of values of the storage
     *  @param[in] fileType   specifies the type of representation (e.g. binary)
     *  @param[in] fileName   is the name of the header file to write
     *  @param[in] size       stands for number of partitions, is CSR storage is part of a distributed matrix
     *  @param[in] rank       stands for the rank of the partition.
     *
     *  Note: for distributed CSR matrices each partition has its own header file.
     *        number of rows and non-zero values are the values for each local part.
     */
    static void writeCSRHeader(
        const IndexType numRows,
        const IndexType numValues,
        const File::FileType& fileType,
        const std::string& fileName,
        const PartitionId size,
        const PartitionId rank );

    static void readCSRHeader(
        IndexType& numRows,
        IndexType& numColumns,
        IndexType& numValues,
        PartitionId& size,
        PartitionId& rank,
        File::FileType& fileType,
        const std::string& fileName );

    /** This method determines some properties of an input file containing matrix data.
     *
     *  param[out] fileType is the type of the file.
     *  param[out] size is the number of files that contain the distributed matrix
     *  param[out] baseFileName is the base name of the file, suffix is given by file type
     *  param[in]  fileName is the name of the input file to read
     */
    static void getFileInfo(
        File::FileType& fileType,
        PartitionId& size,
        std::string& baseFileName,
        const std::string& fileName );

    /** Help routine that determines the availability of a given file by its name. */

    static bool fileExists( const std::string& fileName );

protected:

    /** Logger for this class */

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private    :

    static const int mIversion; //<! unique identification for version
};

/** TODO[doxy] Complete Description.
 *
 * @tparam ValueType is the type of the matrix values.
 */
template<typename ValueType>
class LAMA_DLL_IMPORTEXPORT StorageIO: public _StorageIO
{
public:

    /** @brief General version of writing CSR storage to a file.
     *
     *  @param[in] size             header information in case of distributed matrices
     *  @param[in] rank             header information in case of distributed matrices
     *  @param[in] csrIA            CSR ia array
     *  @param[in] numColumns       number of (global) columns
     *  @param[in] csrJA            CSR ja array
     *  @param[in] csrValues        CSR values array
     *  @param[in] fileName         name of output file
     *  @param[in] fileType         type of file
     *  @param[in] dataType         specifies precision of real values
     *  @param[in] indexDataTypeIA  specifies precision of IA array
     *  @param[in] indexDataTypeJA  specifies precision of JA array
     */
    static void writeCSRToFile(
        const PartitionId size,
        const PartitionId rank,
        const LAMAArray<IndexType>& csrIA,
        const IndexType numColumns,
        const LAMAArray<IndexType>& csrJA,
        const LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const File::FileType& fileType,
        const File::DataType& dataType,
        const File::IndexDataType indexDataTypeIA /*=LONG*/,
        const File::IndexDataType indexDataTypeJA /*=LONG*/
        );

    /** @brief Writing CSR storage to a formatted file.
     *
     *  @param[in] csrIA            CSR ia array
     *  @param[in] csrJA            CSR ja array
     *  @param[in] csrValues        CSR values array
     *  @param[in] fileName         name of output file
     */
    static void writeCSRToFormattedFile(
        const LAMAArray<IndexType>& csrIA,
        const LAMAArray<IndexType>& csrJA,
        const LAMAArray<ValueType>& csrValues,
        const std::string& fileName );

    /** @brief Writing CSR storage to a binary file.
     *
     *  @param[in] csrIA                CSR ia array
     *  @param[in] csrJA                CSR ja array
     *  @param[in] csrValues            CSR values array
     *  @param[in] fileName             name of output file
     *  @param[in] indexDataTypeSizeIA  TODO[doxy] Complete Description.
     *  @param[in] indexDataTypeSizeJA  TODO[doxy] Complete Description.
     *  @param[in] dataTypeSize         TODO[doxy] Complete Description.
     */
    static void writeCSRToBinaryFile(
        const LAMAArray<IndexType>& csrIA,
        const LAMAArray<IndexType>& csrJA,
        const LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const long indexDataTypeSizeIA,
        const long indexDataTypeSizeJA,
        const long dataTypeSize );

    /** @brief Writing CSR storage to an xdr file.
     *
     *  @param[in] csrIA                CSR ia array
     *  @param[in] csrJA                CSR ja array
     *  @param[in] csrValues            CSR values array
     *  @param[in] fileName             name of output file
     *  @param[in] indexDataTypeSizeIA  TODO[doxy] Complete Description.
     *  @param[in] indexDataTypeSizeJA  TODO[doxy] Complete Description.
     *  @param[in] dataTypeSize         TODO[doxy] Complete Description.
     */
    static void writeCSRToXDRFile(
        const LAMAArray<IndexType>& csrIA,
        const LAMAArray<IndexType>& csrJA,
        const LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const long indexDataTypeSizeIA,
        const long indexDataTypeSizeJA,
        const long dataTypeSize );

    /** Writing CSR storage to an mm file.
     *
     *  @param[in] csrIA            CSR ia array
     *  @param[in] numColumns       number of (global) columns
     *  @param[in] csrJA            CSR ja array
     *  @param[in] csrValues        CSR values array
     *  @param[in] fileName         name of output file
     *  @param[in] dataType         specifies precision of real values
     */
    static void writeCSRToMMFile(
        const LAMAArray<IndexType>& csrIA,
        const IndexType numColumns,
        const LAMAArray<IndexType>& csrJA,
        const LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const File::DataType& dataType );

    /** @brief Reading a CSR storage from a file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] numColumns      number of columns
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     */
    static void readCSRFromFile(
        LAMAArray<IndexType>& csrIA,
        IndexType& numColumns,
        LAMAArray<IndexType>& csrJA,
        LAMAArray<ValueType>& csrValues,
        const std::string& fileName );

    /** @brief Reading a CSR storage from a formatted file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     *  @param[in]  numRows         number of rows
     */
    static void readCSRFromFormattedFile(
        LAMAArray<IndexType>& csrIA,
        LAMAArray<IndexType>& csrJA,
        LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const IndexType numRows );

    /** @brief Reading a CSR storage from a binary file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     *  @param[in]  numRows         number of rows
     *
     *  Be careful: no implicit conversions are supported here, so
     *  the binary file must contain data of exact the same type as needed.
     */
    static void readCSRFromBinaryFile(
        LAMAArray<IndexType>& csrIA,
        LAMAArray<IndexType>& csrJA,
        LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const IndexType numRows );

    /** @brief Reading a CSR storage from an xdr file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     *  @param[in]  numRows         number of rows
     */
    static void readCSRFromXDRFile(
        LAMAArray<IndexType>& csrIA,
        LAMAArray<IndexType>& csrJA,
        LAMAArray<ValueType>& csrValues,
        const std::string& fileName,
        const IndexType numRows );

    /** @brief Reading a CSR storage from an mm file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] numColumns      number of columns
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     */
    static void readCSRFromMMFile(
        LAMAArray<IndexType>& csrIA,
        IndexType& numColumns,
        LAMAArray<IndexType>& csrJA,
        LAMAArray<ValueType>& csrValues,
        const std::string& fileName );
};

/* -------------------------------------------------------------------------- */

}
#endif // LAMA_CSRSTORAGE_HPP_
