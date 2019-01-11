/**
 * @file FileIOTest.cpp
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
 * @brief Test of all FileIO classes that have been registered in the FileIO factory
 * @author Thomas Brandes
 * @date 11.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/storage/Storages.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/utilskernel.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

#include <scai/testsupport/uniquePath.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

#include <memory>

using std::unique_ptr;

using namespace scai;
using namespace common;
using namespace lama;

using hmemo::HArray;
using utilskernel::HArrayUtils;

using scai::testsupport::uniquePath;
using scai::testsupport::GlobalTempDir;

using boost::test_tools::per_element;

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

    ValueType zero = 0;
    IndexType bound = 1;

    auto values = utilskernel::sparseRandomHArray<ValueType>( numRows * numColumns, zero, 0.2f, bound );

    SCAI_LOG_INFO( logger, "setDenseData, values = " << values )

    storage.assign( DenseStorage<ValueType>( numRows, numColumns, std::move( values ) ) );
  
    // force diagonal elements

    HArray<IndexType> diagPos;
    HArrayUtils::setSequence<IndexType>( diagPos, 0, 1, numRows );
    HArray<ValueType> diagValues( numRows, ValueType( 0 ) );
    storage.fillCOO( diagPos, diagPos, diagValues, common::BinaryOp::ADD );
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

    HArray<IndexType> csrIA( numRows + 1, ia );
    HArray<IndexType> csrJA( numValues, ja );
    HArray<ValueType> csrValues( numValues, ValueType( 1 ) );

    storage.setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void setSymmetricData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 5;
    const IndexType numColumns = 5;
    const IndexType ia[]     = { 0,             3,   4,        6,        8, 9 };
    const IndexType ja[]     = { 0,   2,   3,   1,   2,   0,   0,   3,   4  };
    const ValueType values[] = { 1.1, 0.1, 0.6, 1.2, 1.3, 0.1, 0.6, 1.4, 1.5 };
    const IndexType numValues = ia[numRows];

    HArray<IndexType> csrIA( numRows + 1, ia );
    HArray<IndexType> csrJA( numValues, ja );
    HArray<ValueType> csrValues( numValues, values );

    storage.setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Unsupported )
{
    typedef SCAI_TEST_TYPE ValueType;

    HArray<ValueType> array( 10, 1 );

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

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStorageFormatted_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("FormattedStorage: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIOFormatted: write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )

        try
        {
            csrStorage.writeToFile( fileName, "", ScalarType::INTERNAL, ScalarType::INDEX_TYPE, FileIO::FORMATTED );
        } 
        catch( common::UnsupportedException& )
        {
            // if storage is unssupported just skip this test here

            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        {
            IndexType numRows;
            IndexType numCols;
            IndexType numValues;

            fileIO->readStorageInfo( numRows, numCols, numValues, fileName );

            SCAI_LOG_DEBUG( logger, "read storage info for file " << fileName 
                                    << " m = " << numRows << " n = " << numCols << ", nnz = " << numValues )

            BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
        }

        BOOST_REQUIRE_EQUAL( csrStorage.getNumRows(), FileIO::getStorageSize( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE_EQUAL( readStorage.getNumColumns(), csrStorage.getNumColumns() );

        RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

        // due to formatted output we might have lost some precision

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
            {
                BOOST_CHECK( common::Math::abs( csrStorage.getValue( i, j ) - readStorage.getValue( i, j ) ) < eps );
            }
        }

        {
            // read just a contiguous block, without first and last row

            const IndexType n = csrStorage.getNumRows() - 2;

            CSRStorage<ValueType> csrBlock;

            fileIO->readStorage( csrBlock, fileName, 1, n );

            for ( IndexType i = 1; i < csrStorage.getNumRows() - 1; ++i )
            {
                for ( IndexType j = 0; j < csrStorage.getNumColumns(); ++j )
                {
                    BOOST_CHECK( common::Math::abs( csrStorage.getValue( i, j ) - csrBlock.getValue( i - 1, j ) ) < eps );
                }
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

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStorageBinary_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("BinaryStorage: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO(binary): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )

        try
        {
            csrStorage.writeToFile( fileName, "", ScalarType::INTERNAL, ScalarType::INDEX_TYPE, FileIO::BINARY );
        } 
        catch( common::UnsupportedException& )
        {
            // if storage is unssupported just skip this test here

            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

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

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStorageNonSquare_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("NonSquareStorage: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO(default): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )

        try
        {
            csrStorage.writeToFile( fileName );
        } 
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

        // The storage read can have less columns

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE( readStorage.getNumColumns() <= csrStorage.getNumColumns() );

        // We verify that the order of the column indexes has not changed

        csrStorage.sortRows();
        readStorage.sortRows();

        BOOST_TEST( hostReadAccess( csrStorage.getIA() ) == hostReadAccess( readStorage.getIA() ), per_element() );
        BOOST_TEST( hostReadAccess( csrStorage.getJA() ) == hostReadAccess( readStorage.getJA() ), per_element() );
        BOOST_TEST( hostReadAccess( csrStorage.getValues() ) == hostReadAccess( readStorage.getValues() ), per_element() );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );

        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SymmetricStorage )
{
    typedef SCAI_TEST_TYPE ValueType;

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
        setSymmetricData( csrStorage );

        BOOST_CHECK( csrStorage.checkSymmetry() );

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStorageSymmetric_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("SymmetricStorage: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO(default): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )

        try
        {
            csrStorage.writeToFile( fileName );
        } 
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

        // The storage read can have less columns

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE( readStorage.getNumColumns() <= csrStorage.getNumColumns() );

        BOOST_CHECK( readStorage.checkSymmetry() );

        // We verify that the order of the column indexes has not changed

        for ( IndexType i = 0; i < csrStorage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < readStorage.getNumColumns(); ++j )
            {
                BOOST_CHECK_MESSAGE( csrStorage.getValue( i, j ) == readStorage.getValue( i, j ),
                                     "matrix[" << i << ", " << j << "] of " << fileName
                                     << ", written " << csrStorage.getValue( i, j )
                                     << ", read " << readStorage.getValue( i, j ) );
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

            const std::string typeName = TypeTraits<ValueType>::id();
            const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStorageEmpty_" + typeName) + fileSuffix;
            BOOST_TEST_MESSAGE("EmptyStorage: fileName = " << fileName);

            m.clear();

            try
            {
                m.writeToFile( fileName );
            } 
            catch( common::UnsupportedException& )
            {
                continue;
            }

            BOOST_CHECK( FileIO::fileExists( fileName ) );

            setNonSquareData( m );  // just set dummy data to see it will be rewritten

            m.readFromFile( fileName );

            // The storage read can have less columns

            BOOST_CHECK_EQUAL( IndexType( 0 ), m.getNumColumns() );
            BOOST_CHECK_EQUAL( IndexType( 0 ), m.getNumRows() );

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

        auto array = utilskernel::randomHArray<ValueType>( N, 100 );

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outArrayFormatted_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("FormattedArray: fileName = " << fileName);

        SCAI_LOG_DEBUG( logger, "FileIO(formatted): write this array: " << array << " via " << *fileIO << " to " << fileName )

        try
        {
            fileIO->writeArray( array, fileName );
        }
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        // Test the method readArrayInfo for
        {
            IndexType size;
            fileIO->readArrayInfo( size, fileName );

            SCAI_LOG_DEBUG( logger, "readArrayInfo " << fileName << ": size = " << size )

            BOOST_CHECK_EQUAL( N, size );
        }

        // Test the method readArray

        RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

        {
            HArray<ValueType> inArray;

            fileIO->readArray( inArray, fileName );

            SCAI_LOG_DEBUG( logger, "Read from file: " << inArray )

            BOOST_REQUIRE_EQUAL( N, inArray.size() );

            RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

            for ( IndexType i = 0; i < N; ++i )
            {
                ValueType expectedVal = array[i];
                ValueType readVal = inArray[i];

                // Due to the formatted output there might be a loss of precision

                BOOST_CHECK( common::Math::abs( expectedVal - readVal ) < eps );
            }
        }

        // Test the method readArray with first, n arguments

        {
            HArray<ValueType> inArray;

            BOOST_CHECK_THROW(
            {
                fileIO->readArray( inArray, fileName, 0, 5 * N );
            }, common::Exception );

            fileIO->readArray( inArray, fileName, 1, N - 2 );

            SCAI_LOG_DEBUG( logger, "Read block from file: " << inArray )

            BOOST_REQUIRE_EQUAL( inArray.size(), N - 2 );

            for ( IndexType i = 1; i < N - 1; ++i )
            {
                ValueType expectedVal = array[i];
                ValueType readVal = inArray[i - 1];

                // Due to the formatted output there might be a loss of precision

                BOOST_CHECK( common::Math::abs( expectedVal - readVal ) < eps );
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

        IndexType range = 1000;

        auto array = utilskernel::randomHArray<ValueType>( N, range );

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outArrayBinary_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("BinaryArray: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO(binary): write this array: " << array << " via " << *fileIO << " to " << fileName )

        try
        {
            fileIO->writeArray( array, fileName );
        }
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        HArray<ValueType> inArray;

        fileIO->readArray( inArray, fileName );

        BOOST_REQUIRE_EQUAL( inArray.size(), array.size() );

        // due to binary output and using same data type there should be no loss

        BOOST_TEST( hostReadAccess( array ) == hostReadAccess( inArray ), per_element() );

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

        HArray<ValueType> array;

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outEmptyArray_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("EmptyArray: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO: write this empty array: " << array << " via " << *fileIO << " to " << fileName )

        try
        {
            fileIO->writeArray( array, fileName );
        }
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        IndexType N = 10;

        HArray<ValueType> inArray( N, ValueType( 1 ) );

        BOOST_CHECK_EQUAL( N, inArray.size() );

        fileIO->readArray( inArray, fileName );

        BOOST_CHECK_EQUAL( IndexType( 0 ), inArray.size() );

        // due to binary output and using same data type there should be no loss

        BOOST_TEST( hostReadAccess( array ) == hostReadAccess( inArray ), per_element() );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( fileName );

        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( writeSparseTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 100;  

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

        float fillRate = 0.05f;   // make it really sparse
        IndexType bound = 1;

        auto array = utilskernel::sparseRandomHArray<ValueType>( N, 0, fillRate, bound );

        HArray<ValueType> sparseArray;
        HArray<IndexType> sparseIndexes;

        ValueType zero = 0;

        utilskernel::HArrayUtils::buildSparseArray( sparseArray, sparseIndexes, array, zero );

        const std::string typeName = TypeTraits<ValueType>::id();
        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outArraySparse_" + typeName) + fileSuffix;
        BOOST_TEST_MESSAGE("writeSparseTest: fileName = " << fileName);

        SCAI_LOG_INFO( logger, "FileIO(sparse): write this sparse array: " << array << " via " << *fileIO << " to " << fileName )

        try
        {
            fileIO->writeSparse( N, sparseIndexes, sparseArray, fileName );
        }
        catch( common::UnsupportedException& )
        {
            continue;
        }

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        // Test the method readArrayInfo for sparse vector data

        {
            IndexType size;
            fileIO->readArrayInfo( size, fileName );
            SCAI_LOG_DEBUG( logger, "readArrayInfo " << fileName << ": size = " << size )
            BOOST_CHECK_EQUAL( N, size );
        }

        // Test the method readSparse

        RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

        {
            IndexType         inN;
            HArray<ValueType> inValues;
            HArray<IndexType> inPos;

            fileIO->readSparse( inN, inPos, inValues, fileName );

            SCAI_LOG_DEBUG( logger, "readSparse( " << fileName << " ): N = " << inN 
                                     << ", indexes = " << inPos << ", vals = " << inValues );

            BOOST_REQUIRE_EQUAL( N, inN );
            BOOST_REQUIRE_EQUAL( sparseArray.size(), inValues.size() );
            BOOST_REQUIRE_EQUAL( sparseIndexes.size(), inPos.size() );

            // indexes must be equal

            BOOST_TEST( hostReadAccess( sparseIndexes ) == hostReadAccess( inPos ), per_element() );

            for ( IndexType i = 0; i < inValues.size(); ++i )
            {
                ValueType expectedVal = sparseArray[i];
                ValueType readVal = inValues[i];

                // Due to the formatted output there might be a loss of precision

                BOOST_CHECK( common::Math::abs( expectedVal - readVal ) < eps );
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

BOOST_AUTO_TEST_CASE( PatternIOTest )
{
    // DataType = Pattern : no I/O of the matrix values

    typedef DefaultReal ValueType;

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

        fileIO->setDataType( common::ScalarType::PATTERN );

        const std::string fileName = uniquePath(GlobalTempDir::getPath(), "outStoragePattern") + fileSuffix;
        BOOST_TEST_MESSAGE("PatternIOTest: fileName = " << fileName);

        CSRStorage<ValueType> csrStorage;
        setNonSquareData( csrStorage );

        SCAI_LOG_INFO( logger, *fileIO << ": write matrix pattern -> " << fileName
                       << ", matrix = " << csrStorage );

        try
        {
            fileIO->writeStorage( csrStorage, fileName );
        }
        catch( common::UnsupportedException& )
        {
            continue;
        }

        CSRStorage<ValueType> readStorage;

        fileIO->readStorage( readStorage, fileName );

        BOOST_REQUIRE_EQUAL( readStorage.getNumRows(), csrStorage.getNumRows() );
        BOOST_REQUIRE_EQUAL( readStorage.getNumValues(), csrStorage.getNumValues() );

        for ( IndexType i = 0; i < readStorage.getNumRows(); ++i )
        {
            // be careful: readStorage.getNumColumns() <= csrStorage.getnumColumns()

            for ( IndexType j = 0; j < readStorage.getNumColumns(); ++j )
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

        fileIO->setDataType( common::ScalarType::INTERNAL );

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
