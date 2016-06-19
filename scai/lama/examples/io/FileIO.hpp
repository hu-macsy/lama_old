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

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/common/ScalarType.hpp>

#include <string>

namespace scai
{

namespace lama
{

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

class FileIO 
{

public:

    /** Write matrix storage into a file.
     *
     *  @param[in] storage is the (local) matrix data that is written
     *  @param[in] fileName is the name of the output file
     *  @param[in] binary if true data is written binary
     *  @param[in] iaType output format used for row indexes
     *  @param[in] jaType output format used for column indexes
     *  @param[in] valuesType output format used for non-zero values
     *
     *  By default all data is written in exactly the same format as it 
     *  is used in the storage. So for a storage with double values 8 bytes 
     *  are needed for each entry, for a storage with float values 4 bytes.
     */

    virtual void writeStorage(
        const _MatrixStorage& storage,
        const std::string& fileName,
        const bool binary = false,
        const common::scalar::ScalarType iaType = common::scalar::INDEX_TYPE,
        const common::scalar::ScalarType jaType = common::scalar::INDEX_TYPE,
        const common::scalar::ScalarType valuesType = common::scalar::INTERNAL ) const = 0;

    /** Write array of arbitrary type into a file.
     *
     *  @param[in] array is the (local) array data that is written
     *  @param[in] fileName is the name of the output file
     *  @param[in] binary if true data is written binary
     *  @param[in] valuesType output format used for values
     */

    virtual void writeArray(
        const hmemo::_HArray& array,
        const std::string& fileName,
        const bool binary = false,
        const common::scalar::ScalarType valuesType = common::scalar::INTERNAL ) const = 0;

    /** Read (local) matrix storage from a file
     *
     *   - implicit conversion by reading from formatted file
     *   - for binary files the type must match the type of storage
     */

    virtual void readStorage(
        _MatrixStorage& storage,
        const std::string& fileName ) const = 0;

    virtual void readArray(
        hmemo::_HArray& array,
        const std::string& fileName ) const = 0;


    /** File suffix can be used to decide about choice of FileIO class */

    virtual std::string getMatrixFileSuffix() const = 0;

    virtual std::string getVectorFileSuffix() const = 0;

    /** Query if formatted or binary IO is supported 
     *
     *  @param[in] binary if true query support for binary, if false query support for formatted
     */

    virtual bool isSupported( const bool binary ) const = 0;

    /** Help routine to delete file and maybe joint files 
     *
     *  @param[in] fileName file to delete
     *  @returns   0 on success
     */
    virtual int deleteFile( const std::string& fileName ) = 0;
};

}

}
