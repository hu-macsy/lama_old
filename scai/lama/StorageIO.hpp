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
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo/HArray.hpp>

#include <scai/lama/io/FileType.hpp>
#include <scai/lama/io/FileStream.hpp>

// std
#include <fstream>

namespace scai
{

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
class COMMON_DLL_IMPORTEXPORT _StorageIO
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

    static void writeMMHeader(
    	const bool& vector,
		const IndexType& numRows,
		const IndexType& numColumns,
		const IndexType& numValues,
		FileStream& fileName,
		const common::scalar::ScalarType& dataType);

    static void readMMHeader(
		IndexType& numRows,
		IndexType& numColumns,
		IndexType& numValues,
		bool& isPattern,
		bool& isSymmetric,
		FileStream& inFile	);


    /** Help routine that determines the availability of a given file by its name. */

    static bool fileExists( const std::string& fileName );

    static bool hasSuffix( const std::string& fileName, const std::string& suffix );
   
protected:

    /** Logger for this class */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    static const int mIversion; //<! unique identification for version
};

/** TODO[doxy] Complete Description.
 *
 * @tparam ValueType is the type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT StorageIO: public _StorageIO
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
     *  @param[in] valuesType       type in which the values should written to the file
     *  @param[in] iaType           type in which the ia array should written to the file
     *  @param[in] jaType           type in which the ja array should written to the file
     *  @param[in] writeBinary      whether the matrix should be written binary
     */
    static void writeCSRToFile(
        const PartitionId size,
        const PartitionId rank,
        const hmemo::HArray<IndexType>& csrIA,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName,
        const File::FileType& fileType,
        const common::scalar::ScalarType& valuesType,
        const common::scalar::ScalarType& iaType,
        const common::scalar::ScalarType& jaType,
        const bool writeBinary = false );

    /** @brief Writing CSR storage to a formatted file.
     *
     *  @param[in] csrIA            CSR ia array
     *  @param[in] csrJA            CSR ja array
     *  @param[in] csrValues        CSR values array
     *  @param[in] fileName         name of output file
     */
    static void writeCSRToFormattedFile(
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName );

    /** @brief Writing CSR storage to a binary file.
     *
     *  @param[in] size                 header information in case of distributed matrices
     *  @param[in] rank                 header information in case of distributed matrices
     *  @param[in] csrIA                CSR ia array
     *  @param[in] csrJA                CSR ja array
     *  @param[in] csrValues            CSR values array
     *  @param[in] fileName             name of output file
     *  @param[in] iaType               Type as which the ja data should be written
     *  @param[in] jaType               Type as which the ia data should be written
     *  @param[in] valuesType           Type as which the value data should be written
     *  @param[in] writeBinary          whether the matrix should be written binary
     */
    static void writeCSRToSAMGFile(
        const PartitionId size,
        const PartitionId rank,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName,
        const common::scalar::ScalarType iaType,
        const common::scalar::ScalarType jaType,
        const common::scalar::ScalarType valuesType,
        const bool writeBinary );

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
        const hmemo::HArray<IndexType>& csrIA,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName,
        const common::scalar::ScalarType& dataType );

    /** @brief Reading a CSR storage from a file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] numColumns      number of columns
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     */
    static void readCSRFromFile(
        hmemo::HArray<IndexType>& csrIA,
        IndexType& numColumns,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName );

    /** @brief Reading a CSR storage from a SAMG file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[out] numColumns      number of columns
     *  @param[in]  fileName        name of input file
     */
    static void readCSRFromSAMGFile(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        IndexType& numColumns,
        const std::string& fileName );

    /** @brief Reading a CSR storage from an mm file.
     *
     *  @param[out] csrIA           CSR ia array
     *  @param[out] numColumns      number of columns
     *  @param[out] csrJA           CSR ja array
     *  @param[out] csrValues       CSR values array
     *  @param[in]  fileName        name of input file
     */
    static void readCSRFromMMFile(
        hmemo::HArray<IndexType>& csrIA,
        IndexType& numColumns,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        const std::string& fileName );


    /** @brief Reading a dense matrix from file
     *
     *  @param[out] data            data
     *  @param[out] numColumns      number of columns
     *  @param[in]  filename        name of input file
     */
    static void readDenseFromFile( hmemo::HArray<ValueType>& data,
                                   IndexType& numColumns,
                                   const std::string& filename );

    /** @brief Writing a dense matrix to a file
     *
     *  @param[in]  data            data
     *  @param[in]  numColumns      number of columns
     *  @param[in]  filename        name of input file
     *  @param[in]  fileType        format of the output file
     *  @param[in]  dataType        type that should be used for writing the file
     *  @param[in]  writeBinary     whether the output file should be write binary or not
     */
    static void writeDenseToFile( const hmemo::HArray<ValueType>& data,
                                  const IndexType& numColumns,
                                  const std::string& filename,
                                  const File::FileType fileType,
                                  const common::scalar::ScalarType dataType,
                                  const bool writeBinary = false );

    /** @brief Reading a dense matrix SAMG file
    *
    *  @param[out] data            data
    *  @param[out] numColumns      number of columns
    *  @param[in]  filename        name of input file
    */
    static void readDenseFromSAMGFile( hmemo::HArray<ValueType>& data,
                                       IndexType& numColumns,
                                       const std::string& filename );

    /** @brief Writing a dense matrix to a SAMG file
     *
     *  @param[in]  data            data
     *  @param[in]  numColumns      number of columns
     *  @param[in]  filename        name of input file
     *  @param[in]  dataType        type that should be used for writing the file
     *  @param[in]  writeBinary     whether the output file should be write binary or not
     */
    static void writeDenseToSAMGFile( const hmemo::HArray<ValueType>& data,
                                      const IndexType& numColumns,
                                      const std::string& filename,
                                      const common::scalar::ScalarType dataType,
                                      const bool writeBinary = false );

    /** @brief Reading a dense matrix from matrix market file
    *
    *  @param[out] data            data
    *  @param[out] numColumns      number of columns
    *  @param[in]  filename        name of input file
    */
    static void readDenseFromMMFile( hmemo::HArray<ValueType>& data,
                                     IndexType& numColumns,
                                     const std::string& filename );

    /** @brief Writing a dense matrix to a matrix market file
     *
     *  @param[in]  data            data
     *  @param[in]  numColumns      number of columns
     *  @param[in]  filename        name of input file
     *  @param[in]  dataType        type that should be used for writing the file
     *  @param[in]  writeBinary     whether the output file should be write binary or not
     */
    static void writeDenseToMMFile( const hmemo::HArray<ValueType>& data,
                                    const IndexType& numColumns,
                                    const std::string& filename,
                                    const common::scalar::ScalarType dataType,
                                    const bool writeBinary = false );
};

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
