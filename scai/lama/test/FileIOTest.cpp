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

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/utilskernel/LArray.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE( FormattedStorage, ValueType, scai_arithmetic_test_types )
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
     
        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outStorageFormatted_" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIOFormatted: write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )
    
        csrStorage.writeToFile( fileName, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::FORMATTED );

        BOOST_CHECK( FileIO::fileExists( fileName ) );

        CSRStorage<ValueType> readStorage;
        readStorage.readFromFile( fileName );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryStorage, ValueType, scai_arithmetic_test_types )
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
     
        std::string typeName = TypeTraits<ValueType>::id();
        std::string fileName = "outStorageBinary" + typeName + fileSuffix;

        SCAI_LOG_INFO( logger, "FileIO(binary): write this storage: " << csrStorage << " via " << *fileIO << " to " << fileName )
    
        csrStorage.writeToFile( fileName, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::BINARY );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( FormattedArray, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryArray, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_SUITE_END();
