/**
 * @file StorageIO.cpp
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
 * @brief Implementation of static IO routines for matrix storage
 * @author Thomas Brandes
 * @date 27.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/StorageIO.hpp>

// others
#include <lama/io/FileIO.hpp>
#include <lama/io/XDRFileStream.hpp>
#include <lama/io/mmio.hpp>

#include <memory/memory.hpp>

#include <lama/exception/Exception.hpp>

#include <lama/openmp/OpenMPCSRUtils.hpp>
#include <tracing/tracing.hpp>

// boost
#include <common/unique_ptr.hpp>
#include <boost/preprocessor.hpp>

using common::unique_ptr;
using common::scoped_array;

using namespace memory;

namespace lama
{

#define VERSION_ID 22

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( _StorageIO::logger, "StorageIO" )

/* -------------------------------------------------------------------------- */

const int _StorageIO::mIversion = 4;

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToFormattedFile(
    const LAMAArray<IndexType>& csrIA,
    const LAMAArray<IndexType>& csrJA,
    const LAMAArray<ValueType>& csrValues,
    const std::string& fileName )
{
    LAMA_REGION( "StorageIO.writeCSRToFormattedFile" )

    IndexType numRows = csrIA.size() - 1;
    IndexType numValues = csrJA.size();

    LAMA_LOG_INFO( logger,
                   "write CSR (#rows = " << numRows << ", #values = " << numValues << ") formatted to file :'" << fileName << "'" )

    //writing matrix data

    std::ofstream amgfile( fileName.c_str(), std::ios::out ); // open .amg

    ReadAccess<IndexType> ia( csrIA );
    ReadAccess<IndexType> ja( csrJA );
    ReadAccess<ValueType> data( csrValues );

    for( IndexType i = 0; i < numRows + 1; ++i )
    {
        IndexType tmp = ia[i] + 1;
        amgfile << tmp << std::endl;
    }

    for( IndexType j = 0; j < numValues; ++j )
    {
        IndexType tmp = ja[j] + 1;
        amgfile << tmp << std::endl;
    }

    for( IndexType j = 0; j < numValues; ++j )
    {
        amgfile << data[j] << std::endl;
    }

    amgfile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromFormattedFile(
    LAMAArray<IndexType>& csrIA,
    LAMAArray<IndexType>& csrJA,
    LAMAArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    LAMA_REGION( "StorageIO.readCSRFromFormattedFile" )

    //Reading matrix data
    std::ifstream amgfile( fileName.c_str(), std::ios::in ); // open .amg

    if( amgfile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    WriteOnlyAccess<IndexType> ia( csrIA, numRows + 1 );

    for( IndexType i = 0; i < numRows + 1; ++i )
    {
        amgfile >> ia[i];
        --ia[i];
    }

    IndexType numValues = ia[numRows];

    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> data( csrValues, numValues );

    for( IndexType j = 0; j < numValues; ++j )
    {
        amgfile >> ja[j];
        --ja[j];
    }

    for( IndexType j = 0; j < numValues; ++j )
    {
        amgfile >> data[j];
    }

    amgfile.close();
}

/* -------------------------------------------------------------------------- */

template<typename FileType,typename DataType,int offset>
static void writeBinaryData( std::fstream& outFile, const DataType data[], const IndexType n )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        outFile.write( reinterpret_cast<const char*>( data ), sizeof(DataType) * n );
        outFile.flush();
        return;
    }

    // allocate buffer for type conversion and/or adding offset

    scoped_array<FileType> buffer( new FileType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] + offset );
    }

    outFile.write( reinterpret_cast<const char*>( buffer.get() ), sizeof(FileType) * n );
    outFile.flush();
    return;
}

/** Help function to determine size of file with CSR data of certain value type
 *
 *  @tparam ValueType type of matrix value, i.e. float or double
 *  @param[in] numRows  number of rows in CSR data
 *  @param[in] numValues number of values in CSR data
 *
 *  @return number of bytes expected for a binary file containing CSR data of this type
 */

template<typename ValueType>
static FileIO::file_size_t expectedCSRFileSize( const IndexType numRows, const IndexType numValues )
{
    FileIO::file_size_t size = sizeof(IndexType) * ( numRows + 1 ); // size of csrIA

    size += sizeof(IndexType) * numValues; // size of csrJA
    size += sizeof(ValueType) * numValues; // size of csrValues

    return size;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromBinaryFile(
    LAMAArray<IndexType>& csrIA,
    LAMAArray<IndexType>& csrJA,
    LAMAArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    LAMA_LOG_INFO( logger,
                   "read CSR<" << common::getScalarType<ValueType>() << "> storage from binary file " << fileName << ", #rows = " << numRows )

    LAMA_REGION( "StorageIO.readCSRFromBinaryFile" )

    FileIO::file_size_t actualSize = FileIO::getFileSize( fileName.c_str() );

    LAMA_LOG_INFO( logger, "CSR binary file " << fileName << " has size = " << actualSize )

    std::fstream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    WriteOnlyAccess<IndexType> ia( csrIA, numRows + 1 );

    // read offset array IA

    FileIO::readBinaryData<IndexType,IndexType, -1>( inFile, ia.get(), numRows + 1 );

    IndexType numValues = ia[numRows];

    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> values( csrValues, numValues );

    FileIO::readBinaryData<IndexType,IndexType, -1>( inFile, ja.get(), numValues );

    // now we can calculate the expected size

    FileIO::file_size_t expectedSize = expectedCSRFileSize<ValueType>( numRows, numValues );

    if( expectedSize == actualSize )
    {
        LAMA_LOG_INFO( logger, "read binary data of type " << csrValues.getValueType() << ", no conversion" )
        FileIO::readBinaryData<ValueType,ValueType,0>( inFile, values.get(), numValues );
    }
    else if( actualSize == expectedCSRFileSize<float>( numRows, numValues ) )
    {
        LAMA_LOG_WARN( logger, "read binary data of type float, conversion to " << csrValues.getValueType() )
        FileIO::readBinaryData<float,ValueType,0>( inFile, values.get(), numValues );
    }
    else if( actualSize == expectedCSRFileSize<double>( numRows, numValues ) )
    {
        LAMA_LOG_WARN( logger, "read binary data of type double, conversion to " << csrValues.getValueType() )
        FileIO::readBinaryData<double,ValueType,0>( inFile, values.get(), numValues );
    }
    else if( actualSize == expectedCSRFileSize<ComplexFloat>( numRows, numValues ) )
    {
        LAMA_LOG_WARN( logger, "read binary data of type double, conversion to " << csrValues.getValueType() )
        FileIO::readBinaryData<ComplexFloat,ValueType,0>( inFile, values.get(), numValues );
    }
    else if( actualSize == expectedCSRFileSize<ComplexDouble>( numRows, numValues ) )
    {
        LAMA_LOG_WARN( logger, "read binary data of type double, conversion to " << csrValues.getValueType() )
        FileIO::readBinaryData<ComplexDouble,ValueType,0>( inFile, values.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION(
                        "File " << fileName << " has unexpected file size = " << actualSize << ", #rows = " << numRows << ", #values = " << numValues << ", expected for float = " << expectedCSRFileSize<float>( numRows, numValues ) << ", or expected for double = " << expectedCSRFileSize<double>( numRows, numValues ) )
    }

    inFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename FileType,typename DataType,int offset>
static void writeData( XDRFileStream& outFile, const DataType data[], const IndexType n )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        outFile.write( data, n );
        return;
    }

    // allocate buffer for type conversion

    FileType* buffer = new FileType[n];

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] + offset );
    }

    outFile.write( buffer, n );
}

/* -------------------------------------------------------------------------- */

template<typename FileType,typename DataType,int offset>
static void readData( XDRFileStream& inFile, DataType data[], const IndexType n )
{
    if( ( offset == 0 ) && ( typeid(FileType) == typeid(DataType) ) )
    {
        // no type conversion needed

        inFile.read( data, n );
        return;
    }

    // allocate buffer for type conversion

    FileType* buffer = new FileType[n];

    inFile.read( buffer, n ); // read data as required file type

    for( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<DataType>( buffer[i] + offset );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToXDRFile(
    const LAMAArray<IndexType>& csrIA,
    const LAMAArray<IndexType>& csrJA,
    const LAMAArray<ValueType>& csrValues,
    const std::string& fileName,
    const long indexDataTypeSizeIA,
    const long indexDataTypeSizeJA,
    const long dataTypeSize )
{
    LAMA_REGION( "StorageIO.writeCSRToXDRFile" )

    IndexType numValues = csrJA.size();
    IndexType numRows = csrIA.size() - 1;

    ReadAccess<IndexType> iaRead( csrIA );
    ReadAccess<IndexType> jaRead( csrJA );
    ReadAccess<ValueType> dataRead( csrValues );

    XDRFileStream outFile( fileName.c_str(), std::ios::out );

    //Write m_ia with m_nnu + 1 elements
    long nnu = 1;
    // todo: += ?!
    nnu = static_cast<long>( numRows );
    //writing m_ia
    outFile.write( &nnu );
    outFile.write( &indexDataTypeSizeIA );

    if( indexDataTypeSizeIA == sizeof(IndexType) )
    {
        writeData<IndexType,IndexType,1>( outFile, iaRead.get(), numRows + 1 );
    }
    else if( indexDataTypeSizeIA == TypeTraits<long>::size )
    {
        writeData<long,IndexType,1>( outFile, iaRead.get(), numRows + 1 );
    }
    else if( indexDataTypeSizeIA == TypeTraits<int>::size )
    {
        writeData<int,IndexType,1>( outFile, iaRead.get(), numRows + 1 );
    }

    outFile.write( &indexDataTypeSizeIA );
    outFile.write( &nnu );
    //writing m_ja with m_nna elements
    long nna = 1;
    nna = static_cast<long>( numValues );
    outFile.write( &nna );
    outFile.write( &indexDataTypeSizeJA );

    if( indexDataTypeSizeJA == sizeof(IndexType) )
    {
        writeData<IndexType,IndexType,0>( outFile, jaRead.get(), numValues );
    }
    else if( indexDataTypeSizeJA == TypeTraits<long>::size )
    {
        writeData<long,IndexType,0>( outFile, jaRead.get(), numValues );
    }
    else if( indexDataTypeSizeJA == TypeTraits<int>::size )
    {
        writeData<int,IndexType,0>( outFile, jaRead.get(), numValues );
    }

    outFile.write( &indexDataTypeSizeJA );
    outFile.write( &numValues );
    //writing m_data
    outFile.write( &nna );
    outFile.write( &dataTypeSize );

    if( dataTypeSize == sizeof(ValueType) )
    {
        writeData<ValueType,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<double>::size )
    {
        writeData<double,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<float>::size )
    {
        writeData<float,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexFloat>::size )
    {
        writeData<ComplexFloat,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexDouble>::size )
    {
        writeData<ComplexDouble,ValueType,0>( outFile, dataRead.get(), numValues );
    }

    outFile.write( &dataTypeSize );
    outFile.write( &numValues );
    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromXDRFile(
    LAMAArray<IndexType>& csrIA,
    LAMAArray<IndexType>& csrJA,
    LAMAArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    LAMA_REGION( "StorageIO.readCSRFromXDRFile" )

    XDRFileStream xdrFile( fileName.c_str(), std::ios::in );
    int indexDataTypeSizeIA;
    int indexDataTypeSizeJA;
    int dataTypeSize;

    if( !xdrFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open XDR matrix file." )
    }

    //Read Index Vector m_ia with m_nnu + 1 elements
    int nnu; // long nnu;
    xdrFile.read( &nnu );

    // check for mismatch in header and XDR matrix file.

    LAMA_ASSERT_EQUAL_ERROR( numRows, (IndexType ) nnu )

    xdrFile.read( &indexDataTypeSizeIA );

    WriteOnlyAccess<IndexType> m_ia( csrIA, numRows + 1 );

    if( sizeof(IndexType) == indexDataTypeSizeIA )
    {
        readData<IndexType,IndexType, -1>( xdrFile, m_ia, numRows + 1 );
    }
    else if( indexDataTypeSizeIA == TypeTraits<int>::size )
    {
        readData<int,IndexType, -1>( xdrFile, m_ia, numRows + 1 );
    }
    else if( indexDataTypeSizeIA == TypeTraits<long>::size )
    {
        readData<long,IndexType, -1>( xdrFile, m_ia, numRows + 1 );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid index data type size in file " + fileName )
    }

    int indexDataTypeSizeIACheck;
    xdrFile.read( &indexDataTypeSizeIACheck );

    LAMA_ASSERT_EQUAL_ERROR( indexDataTypeSizeIA, indexDataTypeSizeIACheck )

    int nnuCheck;
    xdrFile.read( &nnuCheck );

    LAMA_ASSERT_EQUAL_ERROR( nnuCheck, numRows )

    IndexType numValues = m_ia[numRows];

    //Read Index Vector m_ja with m_nna elements
    int nna;
    xdrFile.read( &nna );

    LAMA_ASSERT_EQUAL_ERROR( numValues, (IndexType ) nna );

    xdrFile.read( &indexDataTypeSizeJA );

    WriteOnlyAccess<IndexType> m_ja( csrJA, numValues );

    if( sizeof(IndexType) == indexDataTypeSizeJA )
    {
        readData<IndexType,IndexType,0>( xdrFile, m_ja.get(), numValues );
    }
    else if( indexDataTypeSizeJA == TypeTraits<long>::size )
    {
        readData<long,IndexType,0>( xdrFile, m_ja.get(), numValues );
    }
    else if( indexDataTypeSizeJA == TypeTraits<int>::size )
    {
        readData<int,IndexType,0>( xdrFile, m_ja.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid index data type size in file " + fileName )
    }

    int indexDataTypeSizeJACheck;
    xdrFile.read( &indexDataTypeSizeJACheck );

    LAMA_ASSERT_EQUAL_ERROR( indexDataTypeSizeJA, indexDataTypeSizeJACheck )

    int nnaCheck;
    xdrFile.read( &nnaCheck );

    LAMA_ASSERT_EQUAL_ERROR( nnaCheck, numValues )

    //Read Index Vector m_data with m_nna elements

    xdrFile.read( &nnaCheck );

    LAMA_ASSERT_EQUAL_ERROR( nnaCheck, numValues )

    xdrFile.read( &dataTypeSize );

    WriteOnlyAccess<ValueType> m_data( csrValues, numValues );

    if( dataTypeSize == TypeTraits<double>::size )
    {
        readData<double,ValueType,0>( xdrFile, m_data.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<float>::size )
    {
        readData<float,ValueType,0>( xdrFile, m_data.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexFloat>::size )
    {
        readData<ComplexFloat,ValueType,0>( xdrFile, m_data.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexDouble>::size )
    {
        readData<ComplexDouble,ValueType,0>( xdrFile, m_data.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid data type size in file " + fileName )
    }

    int dataTypeSizeCheck;
    xdrFile.read( &dataTypeSizeCheck );

    LAMA_ASSERT_EQUAL_ERROR( dataTypeSize, dataTypeSizeCheck )

    xdrFile.read( &nnaCheck );

    LAMA_ASSERT_EQUAL_ERROR( nnaCheck, nna )

    xdrFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToBinaryFile(
    const LAMAArray<IndexType>& csrIA,
    const LAMAArray<IndexType>& csrJA,
    const LAMAArray<ValueType>& csrValues,
    const std::string& amgFileName,
    const long indexDataTypeSizeIA,
    const long indexDataTypeSizeJA,
    const long dataTypeSize )
{
    LAMA_REGION( "StorageIO.writeCSRToBinaryFile " )

    ReadAccess<IndexType> iaRead( csrIA );
    ReadAccess<IndexType> jaRead( csrJA );
    ReadAccess<ValueType> dataRead( csrValues );

    IndexType numRows = csrIA.size() - 1;
    IndexType numValues = csrJA.size();

    LAMA_LOG_INFO( logger,
                   "writeCSRToBinaryFile ( " << amgFileName << ")" << ", #rows = " << numRows << ", #values = " << numValues )

    std::fstream outFile( amgFileName.c_str(), std::ios::out | std::ios::binary );

    // write ia, add offset 1

    if( indexDataTypeSizeIA == TypeTraits<int>::size || sizeof(long) == TypeTraits<int>::size )
    {
        writeBinaryData<int,IndexType,1>( outFile, iaRead.get(), numRows + 1 );
    }
    else if( indexDataTypeSizeIA == TypeTraits<long>::size )
    {
        writeBinaryData<long,IndexType,1>( outFile, iaRead.get(), numRows + 1 );
    }
    else
    {
        COMMON_THROWEXCEPTION( "(write unformatted) Unknown index data type size of IA." )
    }

    // write m_ja

    if( indexDataTypeSizeJA == TypeTraits<int>::size || sizeof(long) == TypeTraits<int>::size )
    {
        writeBinaryData<int,IndexType,1>( outFile, jaRead.get(), numValues );
    }
    else if( indexDataTypeSizeJA == TypeTraits<long>::size )
    {
        writeBinaryData<long,IndexType,1>( outFile, jaRead.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "(write unformatted) Unknown index data type size of JA." )
    }

    if( dataTypeSize == TypeTraits<double>::size )
    {
        writeBinaryData<double,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<float>::size )
    {
        writeBinaryData<float,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexFloat>::size )
    {
        writeBinaryData<ComplexFloat,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else if( dataTypeSize == TypeTraits<ComplexDouble>::size )
    {
        writeBinaryData<ComplexDouble,ValueType,0>( outFile, dataRead.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unknown data type size." )
    }

    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToMMFile(
    const LAMAArray<IndexType>& csrIA,
    const IndexType numColumns,
    const LAMAArray<IndexType>& csrJA,
    const LAMAArray<ValueType>& csrValues,
    const std::string& fileName,
    const File::DataType& dataType )
{
    LAMA_REGION( "StorageIO.writeCSRToMMFile" )

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    MM_typecode matcode;
    mm_initialize_typecode( &matcode );
    mm_set_matrix( &matcode );
    mm_set_sparse( &matcode );

    if( dataType == File::DOUBLE || dataType == File::FLOAT || dataType == File::INTERNAL )
    {
        mm_set_real( &matcode );
    }
    else if( dataType == File::COMPLEX )
    {
        mm_set_complex( &matcode );
    }
    else if( dataType == File::INTEGER )
    {
        mm_set_integer( &matcode );
    }
    else if( dataType == File::PATTERN )
    {
        mm_set_pattern( &matcode );
    }
    else
    {
        COMMON_THROWEXCEPTION( "SparseMatrix::writeMatrixToMMFile: " "unknown datatype." << dataType )
    }

    std::FILE* file;

    if( !( file = std::fopen( fileName.c_str(), "w+" ) ) )
    {
        COMMON_THROWEXCEPTION( "SparseMatrix::writeMatrixToMMFile: '" + fileName + "' could not be opened." )
    }

    mm_write_banner( file, matcode );
    mm_write_mtx_crd_size( file, numRows, numColumns, numValues );

    if( std::fclose( file ) != 0 )
    {
        COMMON_THROWEXCEPTION( "SparseMatrix::writeMatrixToMMFile: '" + fileName + "' could not be closed." )
    }

    file = 0;
    std::ofstream ofile;
    ofile.open( fileName.c_str(), std::ios::out | std::ios::app );

    if( ofile.fail() )
    {
        COMMON_THROWEXCEPTION( "SparseMatrix>::writeMatrixToMMFile: '" + fileName + "' could not be reopened." )
    }

    ReadAccess<IndexType> ia( csrIA );
    ReadAccess<IndexType> ja( csrJA );
    ReadAccess<ValueType> data( csrValues );

    for( IndexType ii = 0; ii < numRows; ++ii )
    {
        for( IndexType jj = ia[ii]; jj < ia[ii + 1]; ++jj )
        {
            ofile << ii + 1 << " " << ja[jj] + 1;

            if( dataType != File::PATTERN )
            {
                ofile << " " << data[jj];
            }

            ofile << std::endl;
        }
    }

    ofile.close();
}

/* -------------------------------------------------------------------------- */

// Template struct describes one Value in a Sparse Matrix
// used for the mtx Reader
template<typename ValueType>
struct MatrixValue
{
    IndexType i;
    IndexType j;
    ValueType v;

    MatrixValue( IndexType row, IndexType col, ValueType val )
                    : i( row ), j( col ), v( val )
    {
    }
};

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromMMFile(
    LAMAArray<IndexType>& csrIA,
    IndexType& numColumns,
    LAMAArray<IndexType>& csrJA,
    LAMAArray<ValueType>& csrValues,
    const std::string& fileName )
{
    std::FILE* file;
    file = fopen( fileName.c_str(), "r" );

    if( !file )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    MM_typecode matcode;
    int errorCode = mm_read_banner( file, &matcode );

    if( errorCode != 0 )
    {
        COMMON_THROWEXCEPTION( "Could not process Matrix Market banner. Cause: '" << getErrorString( errorCode ) << "'." );
    }

    bool isPattern = mm_is_pattern( matcode );

    if( mm_is_complex( matcode ) )
    {
        COMMON_THROWEXCEPTION( "Unsupported data type in file '" << fileName << "'." )
    }

    if( !mm_is_matrix( matcode ) )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "' did not contain a matrix." )
    }

    if( !mm_is_sparse( matcode ) )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "' did not contain a sparse matrix." )
    }

    /* symmetric matrices: only lower triangular matrix is stored */
    /* skew matrices: symmetric and all diagonal entries are zero */

    bool symmetric = mm_is_symmetric( matcode ) || mm_is_skew( matcode );

    IndexType numRows = 0;
    IndexType numValues = 0;

    errorCode = mm_read_mtx_crd_size( file, &numRows, &numColumns, &numValues );

    if( errorCode != 0 )
    {
        COMMON_THROWEXCEPTION(
                        "Could not read values from file '" << fileName << "'. Cause: '" << getErrorString( errorCode ) << "'." );
    }

    LAMA_LOG_INFO( logger,
                   "mmx values: #rows = " << numRows << ", #cols = " << numColumns << ", #values = " << numValues )

    if( std::fclose( file ) != 0 )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "' could not be closed." )
    }

    file = 0;
    std::ifstream ifile;
    ifile.open( fileName.c_str(), std::ios::in );

    if( ifile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not reopen file '" << fileName << "'." )
    }

    WriteOnlyAccess<IndexType> ia( csrIA, numRows + 1 );
    // initialize ia;
#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows + 1; i++ )
    {
        ia[i] = 0;
    }

    // First reading in the beginning of the rows
    // then reading in the values and columns of the rows
    //Jump to the beginning of the Values
    char c = '%';

    while( c == '%' )
    {
        ifile >> c;
        ifile.ignore( 1024, '\n' );
    }

    std::vector<MatrixValue<ValueType> > values;

    //Set the right size of the Vector
    if( symmetric )
    {
        values.reserve( numValues * 2 - numRows );
    }
    else
    {
        values.reserve( numValues );
    }

    IndexType lines = numValues;
    //Create Input Vector
    MatrixValue<ValueType> val( 0, 0, 0 );

    for( int l = 0; l < lines && !ifile.eof(); ++l )
    {
        //TODO Read Vector !!!
        // read ia
        ifile >> val.i;
        ifile >> val.j;

        if( !isPattern )
        {
            ifile >> val.v;
        }

        ++ia[val.i];

        // if the matrix is symmetric, the value appears in row 'column' again.
        if( symmetric )
        {
            if( val.j != val.i )
            {
                ++ia[val.j];
                ++numValues;

                MatrixValue<ValueType> val_symmetric( val.j - 1, val.i - 1, val.v );
                values.push_back( val_symmetric );
            }
        }

        ifile.ignore( 256, '\n' );
        val.i--;
        val.j--;
        values.push_back( val );
    }

    if( ifile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    ifile.close();
    ifile.close();
    // Create csrJA, csrValues
    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> data( csrValues, numValues );

    //create absolute Values of ia

    for( int i = 1; i < numRows; i++ )
    {
        ia[i] += ia[i - 1];
        LAMA_LOG_INFO( logger, "offset[" << i << "] = " << ia[i] )
    }

    ia[numRows] = numValues;
    //initialize ia and data
#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < numValues; i++ )
    {
        ja[i] = -1;
        data[i] = 0;
    }

    for( IndexType elem = 0; elem < numValues; elem++ )
    {
        MatrixValue<ValueType> value = values[elem];
        IndexType offset = ia[value.i];
        IndexType pos = 0;

        while( ja[offset + pos] != -1 )
        {
            pos++;
        }

        LAMA_LOG_INFO( logger, "added row " << value.i << ", offset = " << offset + pos << ", j = " << value.j )
        ja[offset + pos] = value.j;
        data[offset + pos] = value.v;
    }

    OpenMPCSRUtils::sortRowElements( ja.get(), data.get(), ia.get(), numRows, true );
    // Note: we do care if the matrix has really all diagonal elements available
    LAMA_LOG_INFO( logger,
                   "construct matrix " << numRows << " x " << numColumns << " from CSR arrays, # non-zeros = " << numValues )
}

/* -------------------------------------------------------------------------- */

bool _StorageIO::fileExists( const std::string& fileName )
{
    // open the file for reading -> it exists

    std::ifstream file( fileName.c_str(), std::ios::in );
    return file.good();
}

/* -------------------------------------------------------------------------- */

#define AMG_FILE_MATRIX_HEADER_SUFFIX "frm"
#define MATRIX_MARKET_FILE_SUFFIX     "mtx"

void _StorageIO::getFileInfo(
    File::FileType& fileType,
    PartitionId& np,
    std::string& baseName,
    const std::string& fileName )
{
    // set default values for non-existing / unknown file

    np = 0;

    fileType = File::FORMATTED;

    baseName = "";

    size_t pos = fileName.find_last_of( "." );

    if( pos == std::string::npos )
    {
        // no suffix available, give it a try with <filename>.amg

        getFileInfo( fileType, np, baseName, fileName + "." + AMG_FILE_MATRIX_HEADER_SUFFIX );

        if( np >= 1 )
        {
            return;
        }

        getFileInfo( fileType, np, baseName, fileName + ".0." + AMG_FILE_MATRIX_HEADER_SUFFIX );

        if( np >= 1 )
        {
            return;
        }

        getFileInfo( fileType, np, baseName, fileName + "." + MATRIX_MARKET_FILE_SUFFIX );

        return;
    }

    if( !fileExists( fileName ) )
    {
        return;
    }

    baseName = fileName.substr( 0, pos - 1 );

    std::string suffix = fileName.substr( pos + 1 );

    LAMA_LOG_DEBUG( logger, "File info of " << fileName << ": base = " << baseName << ", suffix = " << suffix )

    if( suffix == MATRIX_MARKET_FILE_SUFFIX )
    {
        np = 1;
        fileType = File::MATRIX_MARKET;
    }
    else if( suffix == AMG_FILE_MATRIX_HEADER_SUFFIX )
    {
        // read the header file to find out about the number of partitions

        IndexType numRows;
        IndexType numColumns;
        IndexType numValues;
        PartitionId rank;

        LAMA_LOG_DEBUG( logger, "readCSRHeader " << fileName )

        readCSRHeader( numRows, numColumns, numValues, np, rank, fileType, fileName );

        if( np > 1 )
        {
            // ToDo: base name must be cut
        }
    }
}

/* -------------------------------------------------------------------------- */

void _StorageIO::writeCSRHeader(
    const IndexType numRows,
    const IndexType numValues,
    const File::FileType& fileType,
    const std::string& fileName,
    const PartitionId size,
    const PartitionId rank )
{
    char charFileType;

    switch( fileType )
    {
        case File::BINARY:
            charFileType = 'b';
            break;

        case File::FORMATTED:
            charFileType = 'f';
            break;

        case File::XDR:
            charFileType = 'x';
            break;

        default:
            COMMON_THROWEXCEPTION( "Invalid header file." )
    }

    std::fstream outFile( fileName.c_str(), std::ios::out );

    if( !outFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open matrix header file " + fileName + "." )
    }

    outFile << charFileType;
    outFile << " \t";
    outFile << mIversion;
    outFile << "\n";
    outFile << "\t\t";
    outFile << numValues << "\t";
    outFile << numRows << "\t";
    outFile << VERSION_ID << "\t";
    outFile << size << "\t";
    outFile << rank;
    outFile.close();
}

/* -------------------------------------------------------------------------- */

void _StorageIO::readCSRHeader(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    PartitionId& size,
    PartitionId& rank,
    File::FileType& fileType,
    const std::string& frmFileName )
{
    std::ifstream frmFile( frmFileName.c_str(), std::ios::in ); // open *.frm

    if( !frmFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Could not open header file " + frmFileName + "." )
    }

    int iversion;
    char ch = '!';

    numRows = 0;
    numColumns = 0;
    numValues = 0;

    frmFile >> ch >> iversion;

    switch( ch )
    {
        case 'f':
        {
            fileType = File::FORMATTED;
            break;
        }

        case 'b':
        {
            fileType = File::BINARY;
            break;
        }

        case 'x':
        {
            fileType = File::XDR;
            break;
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Invalid file format: " << ch << "." )
        }
    } //switch (ch)

    if( iversion != mIversion )
    {
        COMMON_THROWEXCEPTION( "Invalid file version: " << iversion << ", should be " << mIversion )
    }

    frmFile >> numValues;
    frmFile >> numRows;

    // AMG CSRSparseMatrix is always square

    numColumns = numRows;

    int id;

    frmFile >> id;

    // not really important, may be warning: LAMA_ASSERT_EQUAL_DEBUG( VERSION_ID, id )

    frmFile >> size;
    frmFile >> rank;

    frmFile.close(); // explicitly, otherwise done by destructor
}

/* -------------------------------------------------------------------------- */

static bool hasSuffix( const std::string& name, const char* suffix )
{
    size_t lenSuffix = strlen( suffix );

    if( name.size() >= lenSuffix )
    {
        return name.substr( name.size() - lenSuffix, lenSuffix ) == suffix;
    }

    return false;
}

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToFile(
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
    )
{
    LAMA_REGION( "StorageIO.writeCSRToFile " )

    std::string fileBaseName;

    if( hasSuffix( fileName, ".frm" ) )
    {
        fileBaseName = fileName.substr( 0, fileName.size() - 4 ).c_str();
    }
    else
    {
        fileBaseName = fileName.c_str();
    }

    // for multiple parititions we generate an own id for each partition

    if( size > 1 )
    {
        char rankstr[10];
        sprintf( rankstr, ".%d", rank );
        fileBaseName += rankstr;
    }

    long dataTypeSize = getDataTypeSize<ValueType>( dataType );
    long indexDataTypeSizeIA = getIndexDataTypeSize( indexDataTypeIA );
    long indexDataTypeSizeJA = getIndexDataTypeSize( indexDataTypeJA );

    LAMA_ASSERT_ERROR( indexDataTypeSizeIA > 0, "indexDataTypeIA = " << indexDataTypeIA << " unsupported" )
    LAMA_ASSERT_ERROR( indexDataTypeSizeJA > 0, "indexDataTypeJA = " << indexDataTypeJA << " unsupported" )
    LAMA_ASSERT_ERROR( dataTypeSize >= 0, "dataTypeSize = " << dataTypeSize << " unsupported" )

    switch( fileType )
    {
        case File::FORMATTED:
        {
            writeCSRToFormattedFile( csrIA, csrJA, csrValues, fileBaseName + ".amg" );

            break;
        }

        case File::BINARY:
        {
            writeCSRToBinaryFile( csrIA, csrJA, csrValues, fileBaseName + ".amg", indexDataTypeSizeIA,
                                  indexDataTypeSizeJA, dataTypeSize );
            break;
        }

        case File::XDR:
        {
            writeCSRToXDRFile( csrIA, csrJA, csrValues, fileBaseName + ".amg", indexDataTypeSizeIA, indexDataTypeSizeJA,
                               dataTypeSize );
            break;
        }

        case File::MATRIX_MARKET:
        {
            std::string name = fileName;

            if( fileName.substr( fileName.size() - 4, 4 ) != ".mtx" )
            {
                name += ".mtx";
            }

            writeCSRToMMFile( csrIA, numColumns, csrJA, csrValues, name, dataType );
            return;
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Unknown file type definition." )
        }
    } //switch(fileType)

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    writeCSRHeader( numRows, numValues, fileType, fileBaseName + ".frm", size, rank );
}

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromFile(
    LAMAArray<IndexType>& csrIA,
    IndexType& numColumns,
    LAMAArray<IndexType>& csrJA,
    LAMAArray<ValueType>& csrValues,
    const std::string& fileName )
{
    LAMA_LOG_INFO( logger, "read CSR matrix data from file: '" << fileName << "'." )

    LAMA_REGION( "StorageIO.readCSRFromFile" )

    std::string suffix;

    if( fileName.size() >= 4 )
    {
        suffix = fileName.substr( fileName.size() - 4, 4 );
    }

    if( suffix == ".mtx" )
    {
        readCSRFromMMFile( csrIA, numColumns, csrJA, csrValues, fileName );
        return;
    }

    std::string baseFileName = fileName;

    if( suffix == ".frm" )
    {
        baseFileName = fileName.substr( 0, fileName.size() - 4 );
    }

    File::FileType fileType;

    IndexType numRows;
    IndexType numValues;

    std::string frmFileName = baseFileName + ".frm";
    std::string amgFileName = baseFileName + ".amg";

    PartitionId size = 1;
    PartitionId rank = 0;

    readCSRHeader( numRows, numColumns, numValues, size, rank, fileType, frmFileName );

    LAMA_LOG_INFO( logger,
                   "readCSRHeader( " << frmFileName << " ): " << numRows << " x " << numColumns << ", #values = " << numValues )

    switch( fileType )
    {
        case File::FORMATTED:
        {
            readCSRFromFormattedFile( csrIA, csrJA, csrValues, amgFileName, numRows );
            break;
        }

        case File::BINARY:
        {
            // Attention: no type conversion here, so data sizes must fit
            readCSRFromBinaryFile( csrIA, csrJA, csrValues, amgFileName, numRows );
            break;
        }

        case File::XDR:
        {
            readCSRFromXDRFile( csrIA, csrJA, csrValues, amgFileName, numRows );
            break;
        }

        default:
            COMMON_THROWEXCEPTION( "Read storage file: unknown file type = " << fileType )
    }
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_STORAGE_IO_INSTANTIATE(z, I, _)                              \
    template class COMMON_DLL_IMPORTEXPORT StorageIO<ARITHMETIC_TYPE##I> ;

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_STORAGE_IO_INSTANTIATE, _ )

#undef LAMA_STORAGE_IO_INSTANTIATE

/* -------------------------------------------------------------------------- */

}// namespace LAMA
