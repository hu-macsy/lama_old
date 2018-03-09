/**
 * @file IOStreamTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Test routines for the class IOStream
 * @author Thomas Brandes
 * @date 23.06.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/utilskernel.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/testsupport/uniquePath.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace common;
using namespace lama;

using hmemo::HArray;
using utilskernel::HArrayUtils;

using scai::testsupport::uniquePath;
using scai::testsupport::GlobalTempDir;

using boost::test_tools::per_element;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( IOStreamTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.IOStreamTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( writeFormatted, ValueType, scai_numeric_test_types )
{
    const ValueType values[] = { 1, 2, 3, -1, 2 };

    const IndexType n = sizeof( values ) / sizeof( ValueType );

    HArray<ValueType> data( n, values );

    int precision = 3;

    const auto testFileName = uniquePath(GlobalTempDir::getPath(), "IOStreamTest.writeFormatted");

    IOStream outFile( testFileName, std::ios::out );
    outFile.writeFormatted( data, precision );
    outFile.close();

    HArray<ValueType> data1;

    IOStream inFile( testFileName, std::ios::in );
    inFile.readFormatted( data1, n );
    inFile.close();

    BOOST_CHECK( HArrayUtils::maxDiffNorm( data, data1 ) < 1e-3 );

    int rc = std::remove( testFileName.c_str() );
    BOOST_CHECK_EQUAL( 0, rc );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryTest, ValueType, scai_numeric_test_types )
{
    const IndexType n = 5;
    const IndexType randomRange = 100;

    auto data = utilskernel::randomHArray<ValueType>( n, randomRange );

    // ScalarType type = TypeTraits<ValueType>::stype;

    const auto testFileName = uniquePath(GlobalTempDir::getPath(), "IOStreamTest.BinaryTest");

    IOStream outFile( testFileName, std::ios::out | std::ios::binary );

    ScalarType stype = ScalarType::PATTERN;   // that should throw an exception

    BOOST_CHECK_THROW(
    {
        outFile.writeBinary( data, stype );
    }, common::Exception );

    stype = ScalarType::INTERNAL;  // use the same type as data

    outFile.writeBinary( data, stype );
    outFile.close();

    HArray<ValueType> data1;

    IOStream inFile( testFileName, std::ios::in | std::ios::binary );

    inFile.seekg( 0, std::ios::end );
    size_t realSize = inFile.tellg();
    inFile.seekg( 0, std::ios::beg );

    BOOST_CHECK_EQUAL( n * sizeof( ValueType ), realSize );

    inFile.readBinary( data1, n, stype );
    inFile.close();

    BOOST_TEST( hostReadAccess( data ) == hostReadAccess( data1 ), per_element() );

    int rc = std::remove( testFileName.c_str() );
    BOOST_CHECK_EQUAL( 0, rc );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryConvertTest, ValueType, scai_numeric_test_types )
{
    // This test writes an array binary in a file but uses ScalarRepType that has highest precision

    const IndexType n = 5;

    auto data = utilskernel::randomHArray<ValueType>( n, 1 );

    ScalarType stype = TypeTraits<ScalarRepType>::stype;

    const auto testFileName = uniquePath(GlobalTempDir::getPath(), "BinaryConvertTest");

    IOStream outFile( testFileName, std::ios::out | std::ios::binary );
    outFile.writeBinary( data, stype );
    outFile.close();

    HArray<ValueType> data1;

    IOStream inFile( testFileName, std::ios::in | std::ios::binary );

    inFile.seekg( 0, std::ios::end );
    size_t realSize = inFile.tellg();
    inFile.seekg( 0, std::ios::beg );

    BOOST_CHECK_EQUAL( n * common::typeSize( stype ), realSize );

    inFile.readBinary( data1, n, stype );
    inFile.close();

    // binary should guarantee exactly same values

    BOOST_TEST( hostReadAccess( data ) == hostReadAccess( data1 ), per_element() );

    int rc = std::remove( testFileName.c_str() );
    BOOST_CHECK_EQUAL( 0, rc );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
