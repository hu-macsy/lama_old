/**
 * @file MATIOStream.hpp
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
 * @brief Derived class of IOStream dedicated to MATLAB format.
 * @author Thomas Brandes
 * @date 23.11.2016
 */
#pragma once

// base class
#include <scai/lama/io/IOStream.hpp>

#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace lama
{

/** Extension of class IOStream for MATLAB format.
 *
 *  These are the extensions:
 *
 *    - read a data element frm the file
 *
 */

class COMMON_DLL_IMPORTEXPORT MATIOStream : public IOStream
{
public:

    /** Constructor of a new stream for any mode and opens it */

    MATIOStream( const std::string& filename, std::ios_base::openmode mode );

    /** Read the header of a MAT file */

    void readMATFileHeader( int& version, IOStream::Endian& indicator );

    /** Write the header of a MATFile */

    void writeMATFileHeader();

    /** Read one complete data element from a MATLAB file in a string. 
     *
     *  @param[out] dataElement new allocated char array with the element data
     *  @return     number of bytes allocated for the read element
     *
     *  The read element needs further parsing to read storage/array data.
     */

    uint32_t readDataElement( common::scoped_array<char>& dataElement );

    /** 
     * @param[in] buffer is the string that contains the data element 
     * @param[out] dataType is the type of the element
     * @param[out] nBytes is the number of bytes used for the data element
     * @param[out] wBytes is the total number of written bytes for the element incl header
     * @returns    pointer where the element data is storead
     */
    static const char* readDataElementHeader( uint32_t& dataType, uint32_t& nBytes, uint32_t& wBytes, const char buffer[] );

    uint32_t writeDataElementHeader( const uint32_t dataType, const uint32_t nBytes, bool dryRun );

    enum {
        MAT_INT8       = 1,
        MAT_UINT8      = 2,
        MAT_INT16      = 3,
        MAT_UINT16     = 4,
        MAT_INT32      = 5,
        MAT_UINT32     = 6,
        MAT_FLOAT      = 7,
        MAT_DOUBLE     = 9,
        MAT_LDOUBLE    = 11,
        MAT_INT64      = 12,
        MAT_UINT64     = 13,
        MAT_MATRIX     = 14,
        MAT_COMPRESSED = 15,
        MAT_UTF8       = 16,
        MAT_UTF16      = 17,
        MAT_UTF32      = 18
    };

    typedef enum {
        MAT_CELL_CLASS     = 1,
        MAT_STRUCT_CLASS   = 2,
        MAT_OBJECT_CLASS   = 3,
        MAT_CHAR_CLASS     = 4,
        MAT_SPARSE_CLASS   = 5,
        MAT_DOUBLE_CLASS   = 6,
        MAT_FLOAT_CLASS    = 7,
        MAT_INT8_CLASS     = 8,
        MAT_UINT8_CLASS    = 9,
        MAT_INT16_CLASS    = 10,
        MAT_UINT16_CLASS   = 11,
        MAT_INT32_CLASS    = 12,
        MAT_UINT32_CLASS   = 13,
        MAT_INT64_CLASS    = 14,
        MAT_UINT64_CLASS   = 15
    } MATClass;

    static common::scalar::ScalarType matlabType2ScalarType( uint32_t dataType );

    static uint32_t scalarType2MatlabType( common::scalar::ScalarType dataType );

    static common::scalar::ScalarType class2ScalarType( uint32_t dataType );

    static MATClass scalarType2Class( common::scalar::ScalarType stype );

    template<typename ValueType>
    uint32_t writeData( const ValueType* data, uint32_t size, bool dryRun );

    template<typename ValueType>
    static uint32_t getData( ValueType* data, uint32_t size, const char* buffer );

    inline uint32_t writeString( const char* name, bool dryRun );

    static inline uint32_t getString( char* name, uint32_t nameSize, const char* buffer );

    uint32_t writeSparseHeader(
        const IndexType m,
        const IndexType n,
        const IndexType nnz,
        const uint32_t nBytes,
        bool isComplex,
        bool dryRun );
    
    uint32_t writeDenseHeader(
        const IndexType m,
        const IndexType n,
        const uint32_t nBytes,
        common::scalar::ScalarType stype,
        bool dryRun );

    static uint32_t getMatrixInfo( MATClass& matClass, IndexType dims[2], IndexType& nnz, bool& isComplex, const char* data, bool isCell = false );

private:

    uint32_t writePadding( uint32_t size, bool dryRun );

}; // class MATIOStream

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
uint32_t MATIOStream::writeData( const ValueType* data, uint32_t size, bool dryRun )
{
    common::scalar::ScalarType stype = common::getScalarType<ValueType>();

    uint32_t dataType  = scalarType2MatlabType( stype );
    uint32_t nBytes    = size * sizeof( ValueType );

    uint32_t wBytes    = writeDataElementHeader( dataType, nBytes, dryRun );

    if ( !dryRun )
    {
        std::fstream::write( reinterpret_cast<const char*>( data ), nBytes );
    }

    wBytes += nBytes;
    wBytes += writePadding( nBytes, dryRun );

    SCAI_LOG_INFO( logger, "writeData<" << stype << ">, size = " << size << ", wBytes = " << wBytes << ", dryun = " << dryRun )

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
uint32_t MATIOStream::getData( ValueType* data, uint32_t size, const char* buffer  )
{
    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    SCAI_LOG_INFO( logger, "getData( " << common::TypeTraits<ValueType>::id() << ", size = " << size << " )" )

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, buffer );

    SCAI_ASSERT_EQ_ERROR( nBytes, size * sizeof( ValueType), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( dataType, scalarType2MatlabType( common::TypeTraits<ValueType>::stype ), "type mismatch" )

    ::memcpy( data, dataPtr, nBytes );

    return wBytes;
}

/* --------------------------------------------------------------------------------- */

uint32_t MATIOStream::writeString( const char* name, bool dryRun )
{
    uint32_t size = strlen( name );

    return writeData( name, size, dryRun );
}

/* --------------------------------------------------------------------------------- */

uint32_t MATIOStream::getString( char* name, uint32_t nameSize, const char* buffer  )
{
    uint32_t dataType;
    uint32_t nBytes;
    uint32_t wBytes;

    const char* dataPtr = readDataElementHeader( dataType, nBytes, wBytes, buffer );

    SCAI_ASSERT_EQ_ERROR( dataType, MAT_INT8, "type mismatch" )
    SCAI_ASSERT_LT_ERROR( nBytes, nameSize - 1, "too long string" )

    ::memcpy( name, dataPtr, nBytes );

    name[nBytes] = '\0';   // finalize string

    return wBytes;
}

} /* end namespace lama */

} /* end namespace scai */
