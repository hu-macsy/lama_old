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

#include <scai/common/OpenMP.hpp>
#include <scai/common/ScalarType.hpp>

// hpp
#include <scai/lama/StorageIO.hpp>

// local library
#include <scai/lama/io/FileStream.hpp>
//#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/XDRFileStream.hpp>
#include <scai/lama/io/mmio.hpp>
#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/mepr/IOWrapper.hpp>


// internal scai libraries
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;
using common::unique_ptr;
using common::scoped_array;
using sparsekernel::OpenMPCSRUtils;

namespace lama
{

#define VERSION_ID 22

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( _StorageIO::logger, "StorageIO" )

/* -------------------------------------------------------------------------- */

const int _StorageIO::mIversion = 4;

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToFormattedFile(
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName )
{
    SCAI_REGION( "StorageIO.writeCSRToFormattedFile" )

    SCAI_LOG_INFO( logger, "write CSR (#rows = " << (csrIA.size() - 1) << ", #values = " << csrJA.size()
                           << ") formatted to file :'" << fileName << "'" )

    //writing matrix data
    FileStream amgfile( fileName, std::ios::out );

    // TODO: do we have to use variable output types here?
    amgfile.write(csrIA, 1, common::TypeTraits<IndexType>::stype, "\n");
    amgfile.write(csrJA, 1, common::TypeTraits<IndexType>::stype, "\n");
    amgfile.write(csrValues, 0, common::TypeTraits<ValueType>::stype, "\n");

    amgfile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromFormattedFile(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    SCAI_REGION( "StorageIO.readCSRFromFormattedFile" )

    //Reading matrix data

    FileStream inFile( fileName, std::ios::in );
    inFile.read( csrIA, numRows + 1, -1, common::TypeTraits<IndexType>::stype );
    ReadAccess<IndexType> ia( csrIA );
    IndexType numValues = ia[numRows];

    inFile.read( csrJA, numValues, -1, common::TypeTraits<IndexType>::stype );
    inFile.read( csrValues, numValues, 0, common::TypeTraits<ValueType>::stype );


//    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
//    WriteOnlyAccess<ValueType> data( csrValues, numValues );
//
//    for( IndexType j = 0; j < numValues; ++j )
//    {
////        amgfile >> ja[j];
//        inFile >> ja[j];
//        --ja[j];
//    }
//
//    for( IndexType j = 0; j < numValues; ++j )
//    {
////        amgfile >> data[j];
//        inFile >> data[j];
//    }
//
////    amgfile.close();
    inFile.close();
}

/* -------------------------------------------------------------------------- */

/** Help function to determine size of file with CSR data of certain value type
 *
 *  @tparam ValueType type of matrix value, i.e. float or double
 *  @param[in] numRows  number of rows in CSR data
 *  @param[in] numValues number of values in CSR data
 *
 *  @return number of bytes expected for a binary file containing CSR data of this type
 */

template<typename ValueType>
static IOUtils::file_size_t expectedCSRFileSize( const IndexType numRows, const IndexType numValues )
{
    IOUtils::file_size_t size = sizeof(IndexType) * ( numRows + 1 ); // size of csrIA

    size += sizeof(IndexType) * numValues; // size of csrJA
    size += sizeof(ValueType) * numValues; // size of csrValues

    return size;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromBinaryFile(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "read CSR<" << common::getScalarType<ValueType>() << "> storage from binary file " << fileName << ", #rows = " << numRows )

    SCAI_REGION( "StorageIO.readCSRFromBinaryFile" )

    IOUtils::file_size_t actualSize = IOUtils::getFileSize( fileName.c_str() );

    SCAI_LOG_INFO( logger, "CSR binary file " << fileName << " has size = " << actualSize )

    std::fstream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    WriteOnlyAccess<IndexType> ia( csrIA, numRows + 1 );

    // read offset array IA

    IOUtils::readBinaryData<IndexType,IndexType>( inFile, ia.get(), numRows + 1, -1 );

    IndexType numValues = ia[numRows];

    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> values( csrValues, numValues );

    IOUtils::readBinaryData<IndexType,IndexType>( inFile, ja.get(), numValues, -1 );

    // now we can calculate the expected size

    IOUtils::file_size_t expectedSize = expectedCSRFileSize<ValueType>( numRows, numValues );

    if ( expectedSize == actualSize )
    {
        SCAI_LOG_INFO( logger, "read binary data of type " << csrValues.getValueType() << ", no conversion" )
        IOUtils::readBinaryData<ValueType,ValueType>( inFile, values.get(), numValues );
    }
    else if( mepr::IOWrapper<ValueType, ARITHMETIC_HOST_LIST>::readBinary( actualSize, inFile, values.get(), numValues ))
    {
        SCAI_LOG_WARN( logger, "read binary data of different type, conversion done" )
    }
    else
    {
        COMMON_THROWEXCEPTION(
                        "File " << fileName << " has unexpected file size = " << actualSize << ", #rows = " << numRows << ", #values = " << numValues 
                         << ", expected for float = " << expectedCSRFileSize<float>( numRows, numValues ) 
                         << ", or expected for double = " << expectedCSRFileSize<double>( numRows, numValues ) )
    }

    inFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToXDRFile(
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName,
    const long indexDataTypeSizeIA,
    const long indexDataTypeSizeJA,
    const long dataTypeSize )
{
    SCAI_REGION( "StorageIO.writeCSRToXDRFile" )

    IndexType numValues = csrJA.size();
    IndexType numRows = csrIA.size() - 1;

    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> iaRead( csrIA, host );
    ReadAccess<IndexType> jaRead( csrJA, host );
    ReadAccess<ValueType> dataRead( csrValues, host );

    XDRFileStream outFile( fileName.c_str(), std::ios::out );

    //Write m_ia with m_nnu + 1 elements
    long nnu = 1;
    // todo: += ?!
    nnu = static_cast<long>( numRows );
    //writing m_ia
    outFile.write( &nnu );
    outFile.write( &indexDataTypeSizeIA );

    if( indexDataTypeSizeIA == sizeof( IndexType ) )
    {
        IOUtils::writeXDR<IndexType,IndexType>( outFile, iaRead.get(), numRows + 1, 1 );
    }
    else if( indexDataTypeSizeIA == sizeof( long ) )
    {
        IOUtils::writeXDR<long,IndexType>( outFile, iaRead.get(), numRows + 1, 1 );
    }
    else if( indexDataTypeSizeIA == sizeof ( int ) )
    {
        IOUtils::writeXDR<int,IndexType>( outFile, iaRead.get(), numRows + 1, 1 );
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
        IOUtils::writeXDR<IndexType,IndexType>( outFile, jaRead.get(), numValues );
    }
    else if( indexDataTypeSizeJA == sizeof(long) )
    {
        IOUtils::writeXDR<long,IndexType>( outFile, jaRead.get(), numValues );
    }
    else if( indexDataTypeSizeJA == sizeof(int) )
    {
        IOUtils::writeXDR<int,IndexType>( outFile, jaRead.get(), numValues );
    }

    outFile.write( &indexDataTypeSizeJA );
    outFile.write( &numValues );
    //writing m_data
    outFile.write( &nna );
    outFile.write( &dataTypeSize );

    if ( dataTypeSize == sizeof( ValueType ) )
    {
        IOUtils::writeXDR<ValueType, ValueType>( outFile, dataRead.get(), numValues );
    }
    else if( mepr::IOWrapper<ValueType, ARITHMETIC_HOST_LIST>::writeXDR( dataTypeSize, outFile, dataRead.get(), numValues ) )
    {
        SCAI_LOG_INFO( logger, "writeXDR conversion neeeded" )
    }
    else
    {
        COMMON_THROWEXCEPTION( "no matching type found" )
    }

    outFile.write( &dataTypeSize );
    outFile.write( &numValues );
    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromXDRFile(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName,
    const IndexType numRows )
{
    SCAI_REGION( "StorageIO.readCSRFromXDRFile" )

    XDRFileStream xdrFile( fileName.c_str(), std::ios::in );
    int indexDataTypeSizeIA;
    int indexDataTypeSizeJA;
    int dataTypeSize;

    if( !xdrFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open XDR matrix file." )
    }

    // Read Index Vector m_ia with m_nnu + 1 elements

    int nnu; // long nnu;

    xdrFile.read( &nnu );

    SCAI_ASSERT_EQ_ERROR( numRows, (IndexType ) nnu, "mismatch header and XDR matrix file" )

    xdrFile.read( &indexDataTypeSizeIA );

    WriteOnlyAccess<IndexType> m_ia( csrIA, numRows + 1 );

    if( sizeof(IndexType) == indexDataTypeSizeIA )
    {
        IOUtils::readXDR<IndexType,IndexType>( xdrFile, m_ia, numRows + 1, -1 );
    }
    else if( indexDataTypeSizeIA == sizeof( int ) )
    {
        IOUtils::readXDR<int,IndexType>( xdrFile, m_ia, numRows + 1, -1 );
    }
    else if( indexDataTypeSizeIA == sizeof( long ) )
    {
        IOUtils::readXDR<long,IndexType>( xdrFile, m_ia, numRows + 1, -1 );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid index data type size in file " + fileName )
    }

    int indexDataTypeSizeIACheck;
    xdrFile.read( &indexDataTypeSizeIACheck );

    SCAI_ASSERT_EQ_ERROR( indexDataTypeSizeIA, indexDataTypeSizeIACheck, "size mismatch" )

    int nnuCheck;
    xdrFile.read( &nnuCheck );

    SCAI_ASSERT_EQ_ERROR( nnuCheck, numRows, "")

    IndexType numValues = m_ia[numRows];

    //Read Index Vector m_ja with m_nna elements
    int nna;
    xdrFile.read( &nna );

    SCAI_ASSERT_EQ_ERROR( numValues, (IndexType ) nna, "size mismatch" );

    xdrFile.read( &indexDataTypeSizeJA );

    WriteOnlyAccess<IndexType> m_ja( csrJA, numValues );

    if( sizeof(IndexType) == indexDataTypeSizeJA )
    {
        IOUtils::readXDR<IndexType,IndexType>( xdrFile, m_ja.get(), numValues );
    }
    else if( indexDataTypeSizeJA == sizeof( long ) )
    {
        IOUtils::readXDR<long,IndexType>( xdrFile, m_ja.get(), numValues );
    }
    else if( indexDataTypeSizeJA == sizeof( int ) )
    {
        IOUtils::readXDR<int,IndexType>( xdrFile, m_ja.get(), numValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid index data type size in file " + fileName )
    }

    int indexDataTypeSizeJACheck;
    xdrFile.read( &indexDataTypeSizeJACheck );

    SCAI_ASSERT_EQ_ERROR( indexDataTypeSizeJA, indexDataTypeSizeJACheck, "size mismatch" )

    int nnaCheck;
    xdrFile.read( &nnaCheck );

    SCAI_ASSERT_EQ_ERROR( nnaCheck, numValues, "size mismatch" )

    //Read Index Vector m_data with m_nna elements

    xdrFile.read( &nnaCheck );

    SCAI_ASSERT_EQ_ERROR( nnaCheck, numValues, "size mismatch" )

    xdrFile.read( &dataTypeSize );

    WriteOnlyAccess<ValueType> m_data( csrValues, numValues );

    if( mepr::IOWrapper<ValueType, ARITHMETIC_HOST_LIST>::readXDR( dataTypeSize, xdrFile, m_data.get(), numValues ) )
    {
        SCAI_LOG_TRACE( logger, "read xdr file")
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid data type size in file " + fileName )
    }

    int dataTypeSizeCheck;

    xdrFile.read( &dataTypeSizeCheck );

    SCAI_ASSERT_EQ_ERROR( dataTypeSize, dataTypeSizeCheck, "size mismatch" )

    xdrFile.read( &nnaCheck );

    SCAI_ASSERT_EQ_ERROR( nnaCheck, nna, "size mismatch" )

    xdrFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToBinaryFile(
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& amgFileName,
    const common::scalar::ScalarType iaType,
    const common::scalar::ScalarType jaType,
    const common::scalar::ScalarType valuesType )
{
    SCAI_REGION( "StorageIO.writeCSRToBinaryFile " )

    SCAI_LOG_INFO( logger, "writeCSRToBinaryFile ( " << amgFileName << ")" << ", #rows = " << csrIA.size()-1
                           << ", #values = " << csrJA.size() )

    FileStream outFile( amgFileName, std::ios::out | std::ios::trunc | std::ios::binary );
    outFile.write<IndexType>( csrIA, 1, jaType );
    outFile.write<IndexType>( csrJA, 1, iaType );
    outFile.write<ValueType>( csrValues, 0, valuesType );
    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToMMFile(
    const HArray<IndexType>& csrIA,
    const IndexType numColumns,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName,
    const common::scalar::ScalarType& dataType )
{
    SCAI_REGION( "StorageIO.writeCSRToMMFile" )

    // TODO: different output type possible?

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    writeMMHeader( false, numRows, numColumns, numValues, fileName, dataType );

    FileStream outFile( fileName, std::ios::out | std::ios::app );

    // output code runs only for host context
    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( csrIA, host );
    ReadAccess<IndexType> ja( csrJA, host );
    ReadAccess<ValueType> data( csrValues, host );

    for( IndexType ii = 0; ii < numRows; ++ii )
    {
        for( IndexType jj = ia[ii]; jj < ia[ii + 1]; ++jj )
        {
            outFile << ii + 1 << " " << ja[jj] + 1;

            if( dataType != common::scalar::PATTERN )
            {
                outFile << " " << data[jj];
            }

            outFile << std::endl;
        }
    }

    outFile.close();
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
    HArray<IndexType>& csrIA,
    IndexType& numColumns,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName )
{
    bool isSymmetric, isPattern;
    IndexType numRows, numValues;

    readMMHeader( numRows, numColumns, numValues, isPattern, isSymmetric, fileName );

    std::ifstream ifile;
    ifile.open( fileName.c_str(), std::ios::in );

    if( ifile.fail() )
    {
    	SCAI_LOG_DEBUG( logger, "Could not reopen file " )
        COMMON_THROWEXCEPTION( "Could not reopen file '" << fileName << "'." )
    }

    ContextPtr host = Context::getHostPtr();

    WriteOnlyAccess<IndexType> ia( csrIA, host, numRows + 1 );
    // initialize ia;

	#pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
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
    if( isSymmetric )
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

    SCAI_LOG_DEBUG( logger, "beginning read in" )

    std::string line;

    for( int l = 0; l < lines && !ifile.eof(); ++l )
    {
    	std::getline(ifile, line);
    	std::istringstream reader(line);

        reader >> val.i;
        reader >> val.j;

        if( !isPattern )
        {
            reader >> val.v;
        }
	else
        {
            val.v = ValueType( 1.0 );
        }

        ++ia[val.i];

        // if the matrix is symmetric, the value appears in row 'column' again.
        if( isSymmetric )
        {
            if( val.j != val.i )
            {
                ++ia[val.j];
                ++numValues;

                MatrixValue<ValueType> val_symmetric( val.j - 1, val.i - 1, val.v );
                values.push_back( val_symmetric );
            }
        }

//        not needed anymore, because it will be read line by line
//        ifile.ignore( 256, '\n' );
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
        SCAI_LOG_INFO( logger, "offset[" << i << "] = " << ia[i] )
    }

    ia[numRows] = numValues;
    //initialize ia and data
#pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
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

        SCAI_LOG_INFO( logger, "added row " << value.i << ", offset = " << offset + pos << ", j = " << value.j )
        ja[offset + pos] = value.j;
        data[offset + pos] = value.v;
    }

    OpenMPCSRUtils::sortRowElements( ja.get(), data.get(), ia.get(), numRows, true );
    // Note: we do care if the matrix has really all diagonal elements available
    SCAI_LOG_INFO( logger,
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

bool _StorageIO::hasSuffix( const std::string& fileName, const std::string& suffix)
{
    return fileName.size() >= suffix.size() &&
           fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) == 0;
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

    SCAI_LOG_DEBUG( logger, "File info of " << fileName << ": base = " << baseName << ", suffix = " << suffix )

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

        SCAI_LOG_DEBUG( logger, "readCSRHeader " << fileName )

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

    FileStream outFile( fileName, std::ios::out | std::ios::trunc );
    outFile << charFileType << " \t" << mIversion << "\n";
    outFile << "\t\t" << numValues << "\t" << numRows << "\t" << VERSION_ID << "\t" << size << "\t" << rank;
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
    FileStream frmFile( frmFileName, std::ios::in );

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

    // not really important, may be warning: SCAI_ASSERT_EQ_DEBUG( VERSION_ID, id, "version mismatch" )

    frmFile >> size;
    frmFile >> rank;

    frmFile.close(); // explicitly, otherwise done by destructor
}

void _StorageIO::writeMMHeader(
		const bool& vector,
		const IndexType& numRows,
		const IndexType& numColumns,
		const IndexType& numValues,
		const std::string& fileName,
		const common::scalar::ScalarType& dataType )
{
	MM_typecode matcode;
	mm_initialize_typecode( &matcode );
	mm_set_matrix( &matcode );

	if( vector )
	{
		mm_set_dense( &matcode );
	}
	else
	{
		mm_set_sparse( &matcode );
	}

	if( dataType == common::scalar::DOUBLE || dataType == common::scalar::FLOAT || dataType == common::scalar::INTERNAL )
	{
		mm_set_real( &matcode );
	}
	else if( dataType == common::scalar::COMPLEX || dataType == common::scalar::DOUBLE_COMPLEX || dataType == common::scalar::LONG_DOUBLE_COMPLEX)
	{
		mm_set_complex( &matcode );
	}
	else if( dataType == common::scalar::INDEX_TYPE )
	{
		mm_set_integer( &matcode );
	}
	else if( dataType == common::scalar::PATTERN )
	{
		mm_set_pattern( &matcode );
	}
	else
	{
		COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: " "unknown datatype." << dataType )
	}

	std::FILE* file;

	if( !( file = std::fopen( fileName.c_str(), "w+" ) ) )
	{
		COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: '" + fileName + "' could not be opened." )
	}

	mm_write_banner( file, matcode );

	if( vector )
	{
		SCAI_LOG_DEBUG( logger, "write dense --> " << numRows << "x" << numColumns )
		mm_write_mtx_array_size( file, numRows, numColumns );
	} else
	{
		SCAI_LOG_DEBUG( logger, "write sparse --> " << numRows << "x" << numColumns << " (" << numValues << " values)" )
		mm_write_mtx_crd_size( file, numRows, numColumns, numValues );
	}

	if( std::fclose( file ) != 0 )
	{
		COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: '" + fileName + "' could not be closed." )
	}

	file = 0;
}

void _StorageIO::readMMHeader(
		IndexType& numRows,
		IndexType& numColumns,
		IndexType& numValues,
		bool& isPattern,
		bool& isSymmetric,
		const std::string& fileName )
{
	std::FILE* file;
	file = fopen( fileName.c_str(), "r" );

	if( !file )
	{
		SCAI_LOG_DEBUG( logger, "Could not open file " << fileName )
		COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
	}

	MM_typecode matcode;
	int errorCode = mm_read_banner( file, &matcode );

	if( errorCode != 0 )
	{
		SCAI_LOG_DEBUG( logger, "Could not process Matrix Market banner")
		COMMON_THROWEXCEPTION( "Could not process Matrix Market banner. Cause: '" << getErrorString( errorCode ) << "'." );
	}

	isPattern = mm_is_pattern( matcode );

	if( !mm_is_matrix( matcode ) )
	{
		SCAI_LOG_DEBUG( logger, "file did not contain a matrix" )
		COMMON_THROWEXCEPTION( "'" << fileName << "' did not contain a matrix." )
	}

	if( mm_is_sparse( matcode ) )
	{
		SCAI_LOG_DEBUG( logger, "data is sparse" )
		errorCode = mm_read_mtx_crd_size( file, &numRows, &numColumns, &numValues );
	}
	else if( mm_is_dense( matcode ) )
	{
		SCAI_LOG_DEBUG( logger, "data is dense" )
		errorCode = mm_read_mtx_array_size( file, &numRows, &numColumns );
		numValues = numRows * numColumns;
	}

	if( errorCode != 0 )
	{
		SCAI_LOG_DEBUG( logger, "Could not read values from file")
		COMMON_THROWEXCEPTION(
						"Could not read values from file '" << fileName << "'. Cause: '" << getErrorString( errorCode ) << "'." );
	}

	/* symmetric matrices: only lower triangular matrix is stored */
	/* skew matrices: symmetric and all diagonal entries are zero */

	isSymmetric = mm_is_symmetric( matcode ) || mm_is_skew( matcode );

	SCAI_LOG_INFO( logger,
				   "mmx values: #rows = " << numRows << ", #cols = " << numColumns << ", #values = " << numValues )

	if( std::fclose( file ) != 0 )
	{
		COMMON_THROWEXCEPTION( "'" << fileName << "' could not be closed." )
	}
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToFile(
    const PartitionId size,
    const PartitionId rank,
    const HArray<IndexType>& csrIA,
    const IndexType numColumns,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName,
    const File::FileType& fileType,
    const common::scalar::ScalarType& valuesType,
    const common::scalar::ScalarType& iaType,
    const common::scalar::ScalarType& jaType
    )
{
    SCAI_REGION( "StorageIO.writeCSRToFile " )

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

    long dataTypeSize = getDataTypeSize<ValueType>( valuesType );
    long indexDataTypeSizeIA = getIndexDataTypeSize( iaType );
    long indexDataTypeSizeJA = getIndexDataTypeSize( jaType );

    SCAI_ASSERT_ERROR( indexDataTypeSizeIA > 0, "indexDataTypeIA = " << iaType << " unsupported" )
    SCAI_ASSERT_ERROR( indexDataTypeSizeJA > 0, "indexDataTypeJA = " << jaType << " unsupported" )
    SCAI_ASSERT_ERROR( dataTypeSize >= 0, "dataTypeSize = " << dataTypeSize << " unsupported" )

    switch( fileType )
    {
        case File::FORMATTED:
        {
            writeCSRToFormattedFile( csrIA, csrJA, csrValues, fileBaseName + ".amg" );

            break;
        }

        case File::BINARY:
        {
            writeCSRToBinaryFile( csrIA, csrJA, csrValues, fileBaseName + ".amg", iaType,
                                  jaType, valuesType );
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

            if ( !hasSuffix( name, ".mtx" ) )
            {
                name += ".mtx";
            }

            writeCSRToMMFile( csrIA, numColumns, csrJA, csrValues, name, valuesType );
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
    HArray<IndexType>& csrIA,
    IndexType& numColumns,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName )
{
    SCAI_LOG_INFO( logger, "read CSR matrix data from file: '" << fileName << "'." )

    SCAI_REGION( "StorageIO.readCSRFromFile" )

    std::string suffix;
    std::string baseFileName = fileName;
    std::string amgFileName;
    File::FileType fileType;

    IndexType numRows;
    IndexType numValues;

    if( fileName.size() >= 4 )
    {
        suffix = fileName.substr( fileName.size() - 4, 4 );
    }

    if( suffix == ".frm" )
    {
        baseFileName = fileName.substr( 0, fileName.size() - 4 );
    }

    if( suffix == ".mtx" )
    {
    	fileType = File::MATRIX_MARKET;
    }
    else
    {
		std::string frmFileName = baseFileName + ".frm";

		PartitionId size = 1;
		PartitionId rank = 0;

		readCSRHeader( numRows, numColumns, numValues, size, rank, fileType, frmFileName );

		SCAI_LOG_INFO( logger,
                   "readCSRHeader( " << frmFileName << " ): " << numRows << " x " << numColumns << ", #values = " << numValues )

		amgFileName = baseFileName + ".amg";
    }

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

        case File::MATRIX_MARKET:
        {
        	readCSRFromMMFile( csrIA, numColumns, csrJA, csrValues, fileName );
        	break;
        }

        default:
            COMMON_THROWEXCEPTION( "Read storage file: unknown file type = " << fileType )
    }
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StorageIO, ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
