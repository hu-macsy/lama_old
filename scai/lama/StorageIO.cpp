/**
 * @file StorageIO.cpp
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
 * @brief Implementation of static IO routines for matrix storage
 * @author Thomas Brandes
 * @date 27.07.2012
 */

#include <scai/common/OpenMP.hpp>
#include <scai/common/ScalarType.hpp>

// hpp
#include <scai/lama/StorageIO.hpp>

// local library
#include <scai/lama/io/FileStream.hpp>
//#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/mmio.hpp>
#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/mepr/IOWrapper.hpp>


// internal scai libraries
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/exception/IOException.hpp>

namespace scai
{

using namespace hmemo;
using common::unique_ptr;
using common::scoped_array;
using sparsekernel::OpenMPCSRUtils;

namespace lama
{

/** Constant string for Matrix Market Suffix */

static std::string MM_SUFFIX = ".mtx";

/** SAMG file suffixes */

static std::string SAMG_MAT_HEADER_SUFFIX = ".frm";
static std::string SAMG_MAT_DATA_SUFFIX   = ".amg";
static std::string SAMG_VEC_HEADER_SUFFIX = ".frv";
static std::string SAMG_VEC_DATA_SUFFIX   = ".vec";

/** Help routine to get data file name by header file name 
 *
 *  @param[in] headerFileName is the file name of header file
 *  @return    name of the data file
 *
 *  Note: returns same name if no distinction between header and data 
 */
static std::string getDataFileName( const std::string& headerFileName )
{
    std::string result = headerFileName;

    if ( _StorageIO::hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX) )
    {
        size_t len = SAMG_MAT_HEADER_SUFFIX.length();
        result.replace( result.length() - len, len, SAMG_MAT_DATA_SUFFIX );
    }
    else if ( _StorageIO::hasSuffix( headerFileName, SAMG_VEC_HEADER_SUFFIX ) )
    {
        size_t len = SAMG_VEC_HEADER_SUFFIX.length();
        result.replace( result.length() - len, len, SAMG_VEC_DATA_SUFFIX );
    }

    return result;   // same name if no distinction between header and data
}

#define VERSION_ID 22

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( _StorageIO::logger, "StorageIO" )

/* -------------------------------------------------------------------------- */

const int _StorageIO::mIversion = 4;

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromSAMGFile(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    IndexType& numColumns,
    const std::string& headerFileName )
{
    SCAI_ASSERT_ERROR( hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX ), 
                       "SAMG matrix file name '" << headerFileName << "' illegal, must have suffix " << SAMG_MAT_HEADER_SUFFIX )

    SCAI_REGION( "StorageIO.readCSRFromSAMGFile" )
    // start with reading the header
    FileStream inFile( headerFileName, std::ios::in );
    int iversion;
    char fileType = '!';
    IndexType numValues, numRows, id, size, rank;
    numColumns = 0;
    numValues = 0;
    inFile >> fileType >> iversion;

    if ( fileType != 'f' && fileType != 'b' )
    {
        COMMON_THROWEXCEPTION( "Invalid file type" )
    }

    if ( iversion != mIversion )
    {
        COMMON_THROWEXCEPTION( "Invalid file version: " << iversion << ", should be " << mIversion )
    }

    inFile >> numValues;
    inFile >> numRows;
    // AMG is always square?
    numColumns = numRows;
    inFile >> id;
    inFile >> size;
    inFile >> rank;
    inFile.close(); // explicitly, otherwise done by destructor
    // now read *.amg file in correct mode
    std::ios::openmode flags = std::ios::in;

    if ( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( headerFileName );
    inFile.open( dataFileName, flags );
    //TODO: allow different type to be read?
    inFile.read( csrIA, numRows + 1, IndexType( -1 ), common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrJA, numValues, IndexType( -1 ), common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrValues, numValues, ValueType( 0 ), common::TypeTraits<ValueType>::stype, '\n' );
    inFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToSAMGFile(
    const PartitionId size,
    const PartitionId rank,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& headerFileName,
    const common::scalar::ScalarType iaType,
    const common::scalar::ScalarType jaType,
    const common::scalar::ScalarType valuesType,
    const bool writeBinary )
{
    SCAI_ASSERT_ERROR( hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX ), 
                       "SAMG matrix file name '" << headerFileName << "' illegal, must have suffix " << SAMG_MAT_HEADER_SUFFIX )

    SCAI_REGION( "StorageIO.writeCSRToBinaryFile " ) 

    char fileType = writeBinary ? 'b' : 'f';

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    FileStream outFile( headerFileName, std::ios::out | std::ios::trunc );

    outFile << fileType << " \t" << mIversion << "\n";
    outFile << "\t\t" << numValues << "\t" << numRows << "\t" << VERSION_ID << "\t" << size << "\t" << rank;
    outFile.close();
    SCAI_LOG_INFO( logger, "writeCSRToBinaryFile ( " << headerFileName << ")" << ", #rows = " << csrIA.size() - 1
                   << ", #values = " << csrJA.size() )
    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( writeBinary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( headerFileName );

    outFile.open( dataFileName, flags );
    outFile.write<IndexType>( csrIA, 1, jaType, '\n' );
    outFile.write<IndexType>( csrJA, 1, iaType, '\n' );
    outFile.write<ValueType>( csrValues, 0, valuesType, '\n' );
    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToMMFile(
    const HArray<IndexType>& csrIA,
    const IndexType numColumns,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& mmFileName,
    const common::scalar::ScalarType& dataType )
{
    SCAI_ASSERT_ERROR( hasSuffix( mmFileName, MM_SUFFIX ), 
                       "MatrixMarket file name '" << mmFileName << "' illegal, must have suffix " << MM_SUFFIX )

    SCAI_REGION( "StorageIO.writeCSRToMMFile" )
    // TODO: different output type possible?
    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    FileStream outFile( mmFileName, std::ios::out | std::ios::trunc );

    if ( dataType != common::scalar::INTERNAL )
    {
        writeMMHeader( false, numRows, numColumns, numValues, outFile, dataType );
    }
    else
    {
        writeMMHeader( false, numRows, numColumns, numValues, outFile, common::TypeTraits<ValueType>::stype );
    }

    // output code runs only for host context
    ContextPtr host = Context::getHostPtr();
    ReadAccess<IndexType> ia( csrIA, host );
    ReadAccess<IndexType> ja( csrJA, host );
    ReadAccess<ValueType> data( csrValues, host );

    for ( IndexType ii = 0; ii < numRows; ++ii )
    {
        for ( IndexType jj = ia[ii]; jj < ia[ii + 1]; ++jj )
        {
            outFile << ii + 1 << " " << ja[jj] + 1;

            if ( dataType != common::scalar::PATTERN )
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
    const std::string& mmFileName )
{
    SCAI_ASSERT_ERROR( hasSuffix( mmFileName, MM_SUFFIX ), 
                       "MatrixMarket file name '" << mmFileName << "' illegal, must have suffix " << MM_SUFFIX )

    bool isSymmetric, isPattern;
    IndexType numRows, numValues, numValuesFile;
    FileStream inFile( mmFileName, std::ios::in );
    readMMHeader( numRows, numColumns, numValuesFile, isPattern, isSymmetric, inFile );

    SCAI_LOG_DEBUG( logger, "from header: nrows = " << numRows << ", ncols = " << numColumns << ", nnz = " << numValuesFile 
                            << ", isPattern = " << isPattern << ", isSymmetric = " << isSymmetric )

    ContextPtr host = Context::getHostPtr();
    WriteOnlyAccess<IndexType> ia( csrIA, host, numRows + 1 );
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<IndexType> > setVal;
    setVal.getSupportedContext( host );
    // initialize ia;
    setVal[host]( ia, numRows + 1, 0, utilskernel::reduction::COPY );
    std::vector<MatrixValue<ValueType> > values;

    // Set the right size of the Vector
    if ( isSymmetric )
    {
        values.reserve( numValuesFile * 2 - numRows );
    }
    else
    {
        values.reserve( numValuesFile );
    }

    //Create Input Vector
    MatrixValue<ValueType> val( 0, 0, 0 );
    SCAI_LOG_DEBUG( logger, "beginning read in" )
    std::string line;
    numValues = numValuesFile;

    // TODO: there could be comment or blank lines in the file (allowed by specification?), we should over-read them!
    for ( int l = 0; l < numValuesFile && !inFile.eof(); ++l )
    {
        std::getline( inFile, line );

        SCAI_LOG_TRACE( logger, "Line " << l << " of " << numValuesFile << " : " << line )

        std::istringstream reader( line );
        reader >> val.i;
        reader >> val.j;

        if ( !isPattern )
        {
            reader >> val.v;
        }
        else
        {
            val.v = ValueType( 1.0 );
        }

        ++ia[val.i - 1];

        // if the matrix is symmetric, the value appears in row 'column' again.
        if ( isSymmetric )
        {
            if ( val.j != val.i )
            {
                ++ia[val.j - 1];
                ++numValues;
                MatrixValue<ValueType> val_symmetric( val.j - 1, val.i - 1, val.v );
                values.push_back( val_symmetric );
            }
        }

        val.i--;
        val.j--;

        values.push_back( val );
    }

    if ( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << mmFileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << mmFileName << "': invalid file, contains to many elements." )
    }

    inFile.close();

    SCAI_ASSERT_EQUAL( numValues, (IndexType) values.size(), "size mismatch" )

    // Create csrJA, csrValues
    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> data( csrValues, numValues );
    //create absolute Values of ia
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::scan<IndexType> > scan;
    scan.getSupportedContext( host );
    // convert sizes to offset array
    IndexType ntotal = scan[host]( ia, numRows );
    SCAI_ASSERT_EQUAL( ntotal, numValues, "size mismatch" )
    //initialize ia and data
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<ValueType> > setValValues;
    setValValues.getSupportedContext( host );
    setVal[host]( ja, numValues, -1, utilskernel::reduction::COPY );
    setValValues[host]( data, numValues, 0, utilskernel::reduction::COPY );

    for ( IndexType elem = 0; elem < numValues; elem++ )
    {
        MatrixValue<ValueType> value = values[elem];
        IndexType offset = ia[value.i];
        IndexType pos = 0;

        while ( ja[offset + pos] != -1 )
        {
            pos++;
        }

        SCAI_LOG_TRACE( logger, "added row " << value.i << ", offset = " << offset + pos << ", j = " << value.j )
        ja[offset + pos] = value.j;
        data[offset + pos] = value.v;
    }

    static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;
    sortRowElements.getSupportedContext( host );
    // Note: we do care if the matrix has really all diagonal elements available
    sortRowElements[host]( ja, data, ia, numRows, true );
    SCAI_LOG_INFO( logger, "construct matrix " << numRows << " x " << numColumns << " from CSR arrays, # non-zeros = " << numValues )
}

/* -------------------------------------------------------------------------- */

bool _StorageIO::fileExists( const std::string& fileName )
{
    // open the file for reading -> it exists
    std::ifstream file( fileName.c_str(), std::ios::in );
    return file.good();
}

/* -------------------------------------------------------------------------- */

bool _StorageIO::hasSuffix( const std::string& fileName, const std::string& suffix )
{
    size_t suffixSize = suffix.size();

    // Note: hasSuffix( ".frv", ".frv" ) returns also false, avoids empty names

    return fileName.size() > suffixSize &&
           fileName.compare( fileName.size() - suffixSize, suffixSize, suffix ) == 0;
}

/* -------------------------------------------------------------------------- */

std::string _StorageIO::getSuffix( const std::string& fileName )
{
    size_t pos = fileName.find_last_of( "." );

    if ( pos == std::string::npos )
    {
        return "";
    }

    return fileName.substr( pos );
}

/* -------------------------------------------------------------------------- */

int _StorageIO::removeFile( const std::string& fileName )
{
    int rc = std::remove( fileName.c_str() );
 
    if ( rc )
    {
        return rc;
    }

    // delete data files in case of SAMG header files

    std::string dataFileName = getDataFileName( fileName );

    if ( dataFileName != fileName )
    {
        rc = std::remove( dataFileName.c_str() );
    }

    return rc;
}

/* -------------------------------------------------------------------------- */

void _StorageIO::writeMMHeader(
    const bool& vector,
    const IndexType& numRows,
    const IndexType& numColumns,
    const IndexType& numValues,
    FileStream& outFile,
    const common::scalar::ScalarType& dataType )
{
    outFile << "%%MatrixMarket ";

    if ( vector )
    {
        outFile << "vector array ";
    }
    else
    {
        outFile << "matrix coordinate ";
    }

    switch ( dataType )
    {
        case common::scalar::DOUBLE:
        case common::scalar::FLOAT:
            outFile << "real ";
            break;

        case common::scalar::COMPLEX:
        case common::scalar::DOUBLE_COMPLEX:
        case common::scalar::LONG_DOUBLE_COMPLEX:
            outFile << "complex ";
            break;

        case common::scalar::INDEX_TYPE:
            outFile << "integer ";
            break;

        case common::scalar::PATTERN:
            outFile << "pattern ";
            break;

        default:
            COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: " "unknown datatype." << dataType )
    }

    // TODO: Add support for symmetric matrices
    // currently we can only write non-symmetric
    outFile << "general" << std::endl;

    if ( vector )
    {
        outFile << numRows << " " << numColumns << std::endl;
    }
    else
    {
        outFile << numRows << " " << numColumns << " " << numValues << std::endl;
    }
}

/* -------------------------------------------------------------------------- */

void _StorageIO::readMMHeader(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    bool& isPattern,
    bool& isSymmetric,
    FileStream& inFile )
{
    std::string buffer;
    // read %%MatrixMarket
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "%%matrixmarket" )
    {
        COMMON_THROWEXCEPTION( "Given file is no valid matrix market file, expected file to begin with %%MatrixMarket" )
    }

    // read object type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    // check if object type is valid in general
    if ( buffer != "matrix" && buffer != "vector" )
    {
        COMMON_THROWEXCEPTION( "Object type in the given matrix market file is invalid, should be matrix or vector" )
    }

    bool isVector = false;

    if ( buffer == "vector" )
    {
        isVector = true;
    }

    // read file type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    // checkif file type is valid in general
    if ( buffer != "coordinate" && buffer != "array" )
    {
        COMMON_THROWEXCEPTION( "Format type in the given matrix market file is invalid, should be coordinate or array" )
    }

//    if( buffer == "array" )
//    {
//        // TODO: array => dense data, do we need to check this later?
//        COMMON_THROWEXCEPTION( "File type 'array' (dense data) is currently not supported!" )
//    }
    // read data type
    std::getline( inFile, buffer, ' ' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "real" && buffer != "integer" && buffer != "complex" && buffer != "pattern" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be real, integer, complex or pattern" )
    }

    // TODO: allow to return other value types as well => check if the valid type is used
    if ( buffer == "pattern" )
    {
        isPattern = true;
    }
    else
    {
        isPattern = false;
    }

    // read symmetry
    std::getline( inFile, buffer, '\n' );
    std::transform( buffer.begin(), buffer.end(), buffer.begin(), ::tolower );

    if ( buffer != "general" && buffer != "symmetric" && buffer != "skew-symmetric" && buffer != "hermitian" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be general, symmetric, skew-symmetric or hermitian" )
    }

    if ( buffer == "general" )
    {
        isSymmetric = false;
    }
    else
    {
        if ( buffer == "symmetric" )
        {
            isSymmetric = true;
        }
        else
        {
            // TODO: add support!
            COMMON_THROWEXCEPTION( "Symmetry options 'skew-symmetric' and 'hermitian' are currently not supported!" )
        }
    }

    do
    {
        std::getline( inFile, buffer, '\n' );
    }
    while ( buffer.at( 0 ) == '%' );

    std::stringstream bufferSS( buffer );
    bufferSS >> numRows;
    bufferSS >> numColumns;

    // TODO: vector correct here? should it be dense vs sparse?
    if ( !isVector )
    {
        bufferSS >> numValues;
    }
    else
    {
        numValues = numRows * numColumns;
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
    const common::scalar::ScalarType& jaType,
    const bool writeBinary /* = false */ )
{
    SCAI_REGION( "StorageIO.writeCSRToFile " )

    File::FileType myFileType = fileType;
    std::string    myFileName = fileName;   // might be updated for xxx.<rank>.<suffix>

    if ( myFileType == File::DEFAULT )
    {
        // take decision by suffix
        
        if ( hasSuffix( fileName, SAMG_MAT_HEADER_SUFFIX ) )
        {
            myFileType = File::SAMG_FORMAT;
        }
        else if ( hasSuffix( fileName, MM_SUFFIX ) )
        {
            myFileType = File::MATRIX_MARKET;
        }
        else
        {
            SCAI_THROWEXCEPTION( common::IOException, "Cannot determine file type for writing matrix" )
        }
    }

    // for multiple parititions we generate an own id for each partition

    if ( size > 1 )
    {
        // ToDo outfile.xxx  -> outfile.xxx.1, better is outfile.1.xxx
	std::stringstream ss;
        ss << rank;
	myFileName += ss.str();
//        char rankstr[10];
//        sprintf( rankstr, ".%d", rank );
//        myFileName += rankstr;
    }

    switch ( myFileType )
    {
        case File::SAMG_FORMAT:
            writeCSRToSAMGFile( size, rank, csrIA, csrJA, csrValues, myFileName, iaType, jaType, valuesType, writeBinary );
            break;

        case File::MATRIX_MARKET:
            writeCSRToMMFile( csrIA, numColumns, csrJA, csrValues, myFileName, valuesType );
            return;

        default:
            COMMON_THROWEXCEPTION( "Unknown file type definition." )
    }
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

    File::FileType fileType = File::DEFAULT;   // only chose by suffix

    if ( hasSuffix( fileName, SAMG_MAT_HEADER_SUFFIX ) )
    {
        fileType = File::SAMG_FORMAT;
    }
    else if ( hasSuffix( fileName, MM_SUFFIX ) )
    {
        fileType = File::MATRIX_MARKET;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Cannot determine file type for dense vector/matrix from " << fileName
                               << ", must have suffix " << SAMG_MAT_HEADER_SUFFIX << " or " << MM_SUFFIX )

    }

    switch ( fileType )
    {
        case File::SAMG_FORMAT:
            readCSRFromSAMGFile( csrIA, csrJA, csrValues, numColumns, fileName );
            break;

        case File::MATRIX_MARKET:
            readCSRFromMMFile( csrIA, numColumns, csrJA, csrValues, fileName );
            break;

        default:
            COMMON_THROWEXCEPTION( "Read storage file: unknown file type = " << fileType )
    }
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromFile( 
    HArray<ValueType>& data,
    IndexType& numColumns,
    const std::string& fileName )
{
    File::FileType fileType = File::DEFAULT;  // type chosen by name

    if ( hasSuffix( fileName, SAMG_VEC_HEADER_SUFFIX ) )
    {
        fileType = File::SAMG_FORMAT;
    }
    else if ( hasSuffix( fileName, MM_SUFFIX ) )
    {
        fileType = File::MATRIX_MARKET;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Cannot determine file type for dense vector/matrix from " << fileName 
                               << ", must have suffix " << SAMG_VEC_HEADER_SUFFIX << " or " << MM_SUFFIX )
    }

    switch ( fileType )
    {
        case File::SAMG_FORMAT:
            readDenseFromSAMGFile( data, numColumns, fileName );
            break;

        case File::MATRIX_MARKET:
            readDenseFromMMFile( data, numColumns, fileName );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unknown File Type." );
    }
}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToFile( 
    const HArray<ValueType>& data,
    const IndexType& numColumns,
    const std::string& fileName,
    const File::FileType fileType,
    const common::scalar::ScalarType dataType,
    const bool writeBinary /* = false */ )
{
    File::FileType myFileType = fileType;

    if ( myFileType == File::DEFAULT )
    {
        // take decision by suffix

        if ( hasSuffix( fileName, SAMG_MAT_HEADER_SUFFIX ) || 
             hasSuffix( fileName, SAMG_VEC_HEADER_SUFFIX )  )
        {
            myFileType = File::SAMG_FORMAT;
        }
        else if ( hasSuffix( fileName, MM_SUFFIX ) )
        {
            myFileType = File::MATRIX_MARKET;
        }
        else
        {
            COMMON_THROWEXCEPTION( "Cannot determine file type for writing dense data" )
        }
    }

    switch ( myFileType )
    {
        case File::SAMG_FORMAT:
            writeDenseToSAMGFile( data, numColumns, fileName, dataType, writeBinary );
            break;

        case File::MATRIX_MARKET:
            writeDenseToMMFile( data, numColumns, fileName, dataType, writeBinary );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unknown file type definition." );
    }
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromSAMGFile( 
    HArray<ValueType>& data,
    IndexType& numColumns,
    const std::string& headerFileName )
{
    SCAI_ASSERT_ERROR( hasSuffix( headerFileName, SAMG_VEC_HEADER_SUFFIX ), 
                       "SAMG vector file name '" << headerFileName << "' illegal, must have suffix " << SAMG_VEC_HEADER_SUFFIX )

    SCAI_LOG_INFO( logger, "read DenseVector from file '" << headerFileName << "'" )

    char fileType = ' ';
    int dataTypeSize = 0;
    IndexType numRows;
    // start with reading the *.frv header file
    FileStream inFile( headerFileName, std::ios::in );
    inFile >> fileType;
    inFile >> numRows;
    inFile >> dataTypeSize;
    inFile.close();
    // *.frv files can only store vectors
    numColumns = 1;

    if ( fileType != 'b' && fileType != 'f' )
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << headerFileName 
                               << ", fileType = " << fileType << " illegal" )
    }

    // now read *.vec file in correct mode
    std::ios::openmode flags = std::ios::in;

    if ( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( headerFileName );

    inFile.open( dataFileName, flags );
    common::scalar::ScalarType dataType;

    switch ( dataTypeSize )
    {
        // special cases for handling IndexType, int and long, as these are not properly supported yet
        case 4:
            dataType = common::scalar::FLOAT;
            break;

        case 8:
            dataType = common::scalar::DOUBLE;
            break;

        default:
            dataType = common::scalar::UNKNOWN;
            SCAI_LOG_ERROR( logger, "Encountered invalid type size " << dataTypeSize )
    }

    inFile.read( data, numRows, ValueType( 0 ), dataType, '\n' );
    inFile.close();
}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToSAMGFile( 
    const HArray<ValueType>& data,
    const IndexType& numColumns,
    const std::string& headerFileName,
    const common::scalar::ScalarType dataType,
    const bool writeBinary /* = false */ )
{
    SCAI_ASSERT_ERROR( numColumns == 1, "SAMG format can only store dense vectors" )

    SCAI_ASSERT_ERROR( hasSuffix( headerFileName, SAMG_VEC_HEADER_SUFFIX ), 
                       "SAMG vector file name '" << headerFileName << "' illegal, must have suffix " << SAMG_VEC_HEADER_SUFFIX )

    // start by writing the SAMG Vector header file

    char fileType = writeBinary ? 'b' : 'f';

    // TODO: maybe provide this as function at a central place?

    int typeSize = common::mepr::ScalarTypeHelper<SCAI_ARITHMETIC_ARRAY_HOST_LIST>::sizeOf( dataType );

    if ( typeSize == 0 )
    {
        switch ( dataType )
        {
            case common::scalar::INTERNAL:
                typeSize = sizeof( ValueType );
                break;

            default:
                SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << dataType )
                break;
        }
    }

    FileStream outFile( headerFileName, std::ios::out );
    outFile << fileType << std::endl;
    outFile << data.size() << std::endl;
    outFile << typeSize;
    outFile.close();

    // write data into SAMG vector data file

    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( writeBinary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( headerFileName );

    outFile.open( dataFileName, flags );
    outFile.write<ValueType>( data, 0, dataType, '\n' );
    outFile.close();
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromMMFile(
    HArray<ValueType>& data,
    IndexType& numColumns,
    const std::string& mmFileName )

{
    SCAI_ASSERT_ERROR( hasSuffix( mmFileName, MM_SUFFIX ), 
                       "MatrixMarket file name '" << mmFileName << "' illegal, must have suffix " << MM_SUFFIX )

    bool isSymmetric, isPattern;
    IndexType numRows, numValues, i;
    ValueType val;
    std::string line;
    FileStream inFile( mmFileName, std::ios::in );
    readMMHeader( numRows, numColumns, numValues, isPattern, isSymmetric, inFile );
    WriteOnlyAccess<ValueType> vector( data, numValues );
    ValueType* vPtr = vector.get();

    for ( int l = 0; l < numValues && !inFile.eof(); ++l )
    {
        std::getline( inFile, line );
        std::istringstream reader( line );

        if ( isPattern )
        {
            reader >> i;
            val = 1.0;
            i--;
        }
        else
        {
            reader >> val;
            i = l;
        }

        vPtr[i] = val;
    }

    if ( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << mmFileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline( inFile, line );

    if ( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << mmFileName << "': invalid file, contains to many elements." )
    }

    inFile.close();
    SCAI_LOG_INFO( logger, "construct vector " << numRows )
}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToMMFile( const HArray<ValueType>& data,
        const IndexType& numColumns,
        const std::string& mmFileName,
        const common::scalar::ScalarType dataType,
        const bool writeBinary /* = false */ )
{
    SCAI_ASSERT_ERROR( hasSuffix( mmFileName, MM_SUFFIX ), 
                       "MatrixMarket file name '" << mmFileName << "' illegal, must have suffix " << MM_SUFFIX )

    SCAI_ASSERT_ERROR( ! writeBinary, "Matrix market format can not be written binary" );

    FileStream outFile( mmFileName, std::ios::out | std::ios::trunc );

    if ( dataType != common::scalar::INTERNAL )
    {
        writeMMHeader( true, data.size(), numColumns, -1, outFile, dataType );
    }
    else
    {
        writeMMHeader( true, data.size(), numColumns, -1, outFile, common::TypeTraits<ValueType>::stype );
    }

    // output code runs only for host context
    ContextPtr host = Context::getHostPtr();
    ReadAccess<ValueType> dataRead( data, host );

    for ( IndexType i = 0; i < data.size(); ++i )
    {
        outFile << dataRead[i] << std::endl;
    }

    outFile.close();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StorageIO, SCAI_ARITHMETIC_HOST )

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
