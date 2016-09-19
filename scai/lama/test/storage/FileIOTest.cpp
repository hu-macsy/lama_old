/**
 * @file FileIOTest.cpp
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
 * @brief Test of all FileIO classes that have been registered in the FileIO factory
 * @author Thomas Brandes
 * @date 11.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/storage/Storages.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace common;
using namespace lama;
using namespace hmemo;

using utilskernel::LArray;

/** Output files should be deleted unless for debugging it might be useful to check them. */

#define DELETE_OUTPUT_FILES

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( FileIOTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.FileIOTest" );

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 7;
    const IndexType numColumns = 7;

    // values: take numRows x numColums random numbers of required type

    LArray<ValueType> values;
    float fillRate = 0.2;    
    values.setRandom( numRows * numColumns, fillRate );

    ValueType eps = static_cast<ValueType>( 1E-5 );

    // Note: diagonal property of sparse matrices will be set due to square matrix

    storage.setDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void setNonSquareData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4; 
    const IndexType numColumns = 10;
    const IndexType ia[] = { 0,    2,       5, 6,    8 };
    const IndexType ja[] = { 1, 2, 3, 2, 4, 5, 7, 4 };
    const IndexType numValues = ia[numRows];

    LArray<IndexType> csrIA( numRows + 1, ia );
    LArray<IndexType> csrJA( numValues, ja ); 
    LArray<ValueType> csrValues( numValues, ValueType( 1 ) );

    storage.setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Unsupported )
{
    typedef SCAI_TEST_TYPE ValueType;

    LArray<ValueType> array( 10, 1 );

    BOOST_CHECK_THROW(
    {
        FileIO::write( array, "myArray.unsupported_suffix" );
    }, common::Exception );

    BOOST_CHECK_THROW(
    {
        FileIO::read( array, "myArray.unsupported_suffix" );
    }, common::Exception );

    int rc = FileIO::removeFile( "myArray.unsupported_suffix" );

    BOOST_CHECK( rc != 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Settings )
{
    common::Settings::putEnvironment( "SCAI_IO_BINARY", "1" );
    common::Settings::putEnvironment( "SCAI_IO_TYPE_DATA", "_Internal" );

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        common::Settings::putEnvironment( "SCAI_IO_TYPE_INDEX", "_Internal" );

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );

        common::Settings::putEnvironment( "SCAI_IO_TYPE_INDEX", "float2" );

        BOOST_CHECK_THROW (
        {
            unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
        }, common::Exception );
    }

    // remove entries

    unsetenv( "SCAI_IO_BINARY" );
    unsetenv( "SCAI_IO_TYPE_DATA" );
    unsetenv( "SCAI_IO_TYPE_INDEX" );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( WriteAtTest )
{
    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {   
        const std::string& fileSuffix = supportedSuffixes[i];
        
        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );

        SCAI_LOG_INFO( logger, "Suffix " << fileSuffix << ": " << *fileIO )

        std::ostringstream out1;
        out1 << *fileIO;

        BOOST_CHECK( out1.str().length() > 0 );

        std::ostringstream out2;
        fileIO->FileIO::writeAt( out2 );
        BOOST_CHECK( out2.str().length() > 0 );

        // make sure that derived FileIO class has overridden writeAt

        BOOST_CHECK( out1.str() != out2.str() );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FormattedStorage, ValueType, scai_numeric_test_types )
{
    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( !fileIO->isSupportedMode( FileIO::FORMATTED ) )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped, does not support FORMATTED mode" )
            continue;
        }

        if ( fileSuffix != fileIO->getMatrixFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for matrix, is not default matrix suffix" )
            continue;   
        }

        CSRStorage<ValueType> csrStorage;
        setDenseData( csrStorage );
     
        BOOST_CHECK( csrStorage.hasDiagonalProperty() );

        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outStorageFormatted_" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIOFormatted: write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )
    
        csrStorage.writeToFile( fileName, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::FORMATTED );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

        BOOST_CHECK( readStorage.hasDiagonalProperty() );

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE_EQUAL( readStorage.getNumColumns(), csrStorage.getNumColumns() );

        // due to formatted output we might have lost some precision

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
            {
                SCAI_CHECK_CLOSE( csrStorage.getValue( i, j ),
                                  readStorage.getValue( i, j ), 0.01f );
            }
        }

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryStorage, ValueType, scai_numeric_test_types )
{
    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( !fileIO->isSupportedMode( FileIO::BINARY ) )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped, does not support BINARY mode" )
            continue;
        }

        if ( fileSuffix != fileIO->getMatrixFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for matrix, is not default matrix suffix" )
            continue;   
        }

        CSRStorage<ValueType> csrStorage;
        setDenseData( csrStorage );
     
        BOOST_CHECK( csrStorage.hasDiagonalProperty() );

        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outStorageBinary" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO(binary): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )
    
        csrStorage.writeToFile( fileName, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::BINARY );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );
        BOOST_CHECK( readStorage.hasDiagonalProperty() );

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE_EQUAL( readStorage.getNumColumns(), csrStorage.getNumColumns() );

        // due to binary output and using same data type there should be no loss

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
            {
                BOOST_CHECK_EQUAL( csrStorage.getValue( i, j ), readStorage.getValue( i, j ) );
            }
        }

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( NonSquareStorage )
{
    typedef SCAI_TEST_TYPE ValueType;

    // read / write of matrix storage with size 2 x 8
    // Note: for partitioned IO the matrices might be no more square

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( fileSuffix != fileIO->getMatrixFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for matrix, is not default matrix suffix" )
            continue;   
        }

        CSRStorage<ValueType> csrStorage;
        setNonSquareData( csrStorage );

        LArray<IndexType> firstColIndexes1;
        csrStorage.getFirstColumnIndexes( firstColIndexes1 );
     
        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outStorageNonSquare" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO(default): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )
    
        csrStorage.writeToFile( fileName );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

        // The storage read can have less columns 

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE( readStorage.getNumColumns() <= csrStorage.getNumColumns() );

        // We verify that the order of the column indexes has not changed

        LArray<IndexType> firstColIndexes2;
        readStorage.getFirstColumnIndexes( firstColIndexes2 );

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {
            BOOST_CHECK_EQUAL( firstColIndexes1[i], firstColIndexes2[i] );

            for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
            {
                BOOST_CHECK_EQUAL( csrStorage.getValue( i, j ), readStorage.getValue( i, j ) );
            }
        }

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( EmptyStorage )
{
    typedef SCAI_TEST_TYPE ValueType;

    // read / write of matrix storage with size 0 x 0

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( fileSuffix != fileIO->getMatrixFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for matrix, is not default matrix suffix" )
            continue;   
        }

        TypedStorages<ValueType> storages;

        for ( size_t j = 0; j < storages.size(); ++j )
        {
            MatrixStorage<ValueType>& m = *storages[j];

            std::string typeName = TypeTraits<ValueType>::id();
            std::string fileName = "outEmptyStorage" + typeName + fileSuffix;

            m.clear();
            m.writeToFile( fileName );

            BOOST_CHECK( FileIO::fileExists( fileName ) );

            setNonSquareData( m );  // just set dummy data to see it will be rewritten

            m.readFromFile( fileName );

            // The storage read can have less columns 

            BOOST_CHECK_EQUAL( 0, m.getNumColumns() );
            BOOST_CHECK_EQUAL( 0, m.getNumRows() );

#ifdef DELETE_OUTPUT_FILES
            int rc = FileIO::removeFile( fileName );
    
            BOOST_CHECK_EQUAL( rc, 0 );
            BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FormattedArray, ValueType, scai_numeric_test_types )
{
    const IndexType N = 20;

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( !fileIO->isSupportedMode( FileIO::FORMATTED ) )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped, does not support FORMATTED mode" )
            continue;
        }

        if ( fileSuffix != fileIO->getVectorFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for vector, " << fileSuffix 
                                   << " is not default vector suffix" << fileIO->getVectorFileSuffix() )
            continue;   
        }

        fileIO->setMode( FileIO::FORMATTED );

        LArray<ValueType> array;

        array.setRandom( N );
     
        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outArrayFormatted" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO(formatted): write this array: " << array << " via " << *fileIO << " to " << fileName )
    
        fileIO->writeArray( array, fileName );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        LArray<ValueType> inArray;

        fileIO->readArray( inArray, fileName );

        BOOST_REQUIRE_EQUAL( inArray.size(), array.size() );

        // due to binary output and using same data type there should be no loss

        for ( IndexType i = 0; i < N; ++i )
        {
            ValueType expectedVal = array[i];
            ValueType readVal = inArray[i];
            SCAI_CHECK_CLOSE( expectedVal, readVal, 0.01f );
        }


#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryArray, ValueType, scai_numeric_test_types )
{
    const IndexType N = 20;

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( !fileIO->isSupportedMode( FileIO::BINARY ) )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped, does not support BINARY mode" )
            continue;
        }

        if ( fileSuffix != fileIO->getVectorFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for vector, " << fileSuffix 
                                   << " is not default vector suffix" << fileIO->getVectorFileSuffix() )
            continue;   
        }

        fileIO->setMode( FileIO::BINARY );

        LArray<ValueType> array;

        array.setRandom( N );
     
        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outArrayBinary_" + typeName  + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO(binary): write this array: " << array << " via " << *fileIO << " to " << fileName )
    
        fileIO->writeArray( array, fileName );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        LArray<ValueType> inArray;

        fileIO->readArray( inArray, fileName );

        BOOST_REQUIRE_EQUAL( inArray.size(), array.size() );

        // due to binary output and using same data type there should be no loss

        BOOST_CHECK_EQUAL( 0, array.maxDiffNorm( inArray ) );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( EmptyArray, ValueType, scai_numeric_test_types )
{
    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( fileSuffix != fileIO->getVectorFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for vector, " << fileSuffix 
                                   << " is not default vector suffix" << fileIO->getVectorFileSuffix() )
            continue;   
        }

        LArray<ValueType> array;

        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outEmptyArray_" + typeName  + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO: write this empty array: " << array << " via " << *fileIO << " to " << fileName )
    
        fileIO->writeArray( array, fileName );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        IndexType N = 10;

        LArray<ValueType> inArray( N, ValueType( 1 ) );

        BOOST_CHECK_EQUAL( N, inArray.size() );

        fileIO->readArray( inArray, fileName );

        BOOST_CHECK_EQUAL( 0, inArray.size() );

        // due to binary output and using same data type there should be no loss

        BOOST_CHECK_EQUAL( 0, array.maxDiffNorm( inArray ) );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( PatternIOTest )
{
    // DataType = Pattern : no I/O of the matrix values

    typedef RealType ValueType;

    std::vector<std::string> supportedSuffixes;

    FileIO::getCreateValues( supportedSuffixes );

    // loop over all supported suffixes, got them from FileIO factory

    for ( size_t i = 0; i < supportedSuffixes.size(); ++i )
    {
        const std::string& fileSuffix = supportedSuffixes[i];

        unique_ptr<FileIO> fileIO( FileIO::create( fileSuffix ) );
 
        if ( fileSuffix != fileIO->getMatrixFileSuffix() )
        {
            SCAI_LOG_INFO( logger, *fileIO << " skipped for matrix, is not default matrix suffix" )
            continue;   
        }

        fileIO->setDataType( common::scalar::PATTERN );

        std::string fileName = "outStoragePattern" + fileSuffix;

        CSRStorage<ValueType> csrStorage;
        setNonSquareData( csrStorage );

        SCAI_LOG_INFO( logger, *fileIO << ": write matrix pattern -> " << fileName 
                               << ", matrix = " << csrStorage );

        fileIO->writeStorage( csrStorage, fileName );

        CSRStorage<ValueType> readStorage;

        fileIO->readStorage( readStorage, fileName );

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE_EQUAL( readStorage.getNumValues(), csrStorage.getNumValues() );

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {   
            for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
            {   
                if ( csrStorage.getValue( i, j ) != 0 )
                {
                    BOOST_CHECK( readStorage.getValue( i, j ) == ValueType( 1 ) );
                }
                else
                {
                    BOOST_CHECK( readStorage.getValue( i, j ) == 0 );
                }
            }
        }

        // If we want to read a full matrix, the read operation must throw an exception

        fileIO->setDataType( common::scalar::INTERNAL );

        SCAI_LOG_INFO( logger, *fileIO << ": read matrix pattern -> " << fileName 
                               << ", matrix = " << csrStorage );

        BOOST_CHECK_THROW(
        {
            fileIO->readStorage( readStorage, fileName );
        }, common::Exception );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );
    
        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
