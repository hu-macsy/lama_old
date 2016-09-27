/**
 * @file IOStreamTest.cpp
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
 * @brief Test routines for the class IOStream
 * @author Thomas Brandes
 * @date 23.06.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace common;
using namespace lama;
using namespace utilskernel;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( IOStreamTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.IOStreamTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( writeFormatted, ValueType, scai_numeric_test_types )
{
    const IndexType n = 5;

    LArray<ValueType> data( n );

    data[0] = 1.0;
    data[1] = 2.0;
    data[2] = 3.0;
    data[3] = -1.0;
    data[4] = 1;

    int precision = 3;

    IOStream outFile( "tmp.data", std::ios::out );
    outFile.writeFormatted( data, precision );
    outFile.close();

    LArray<ValueType> data1;

    IOStream inFile( "tmp.data", std::ios::in );
    inFile.readFormatted( data1, n );
    inFile.close();

    BOOST_REQUIRE_EQUAL( data.size(), data1.size() );

    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    AbsType diff = common::Math::real( data.maxDiffNorm( data1 ) );

    BOOST_CHECK( diff < 1e-3 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryTest, ValueType, scai_numeric_test_types )
{
    const IndexType n = 5;

    float fillRate = 1.0f;

    LArray<ValueType> data;
    HArrayUtils::setRandom( data, n, fillRate );

    // scalar::ScalarType type = TypeTraits<ValueType>::stype;

    IOStream outFile( "tmp.data", std::ios::out | std::ios::binary );

    scalar::ScalarType stype = scalar::PATTERN;   // that should throw an exception

    BOOST_CHECK_THROW(
    {
        outFile.writeBinary( data, stype );
    }, common::Exception );

    stype = scalar::INTERNAL;  // use the same type as data

    outFile.writeBinary( data, stype );
    outFile.close();

    LArray<ValueType> data1;

    IOStream inFile( "tmp.data", std::ios::in | std::ios::binary );

    inFile.seekg( 0, std::ios::end );
    size_t realSize = inFile.tellg();
    inFile.seekg( 0, std::ios::beg );

    BOOST_CHECK_EQUAL( n * sizeof( ValueType ), realSize );

    inFile.readBinary( data1, n, stype );
    inFile.close();

    BOOST_REQUIRE_EQUAL( data.size(), data1.size() );

    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    AbsType diff = common::Math::real( data.maxDiffNorm( data1 ) );

    BOOST_CHECK_EQUAL( diff, 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( BinaryConvertTest, ValueType, scai_numeric_test_types )
{
    // This test writes an array binary in a file but uses ScalarRepType that has highest precision

    const IndexType n = 5;

    float fillRate = 1.0f;

    LArray<ValueType> data;
    HArrayUtils::setRandom( data, n, fillRate );

    scalar::ScalarType stype = TypeTraits<ScalarRepType>::stype;

    IOStream outFile( "tmp.data", std::ios::out | std::ios::binary );
    outFile.writeBinary( data, stype );
    outFile.close();

    LArray<ValueType> data1;

    IOStream inFile( "tmp.data", std::ios::in | std::ios::binary );

    inFile.seekg( 0, std::ios::end );
    size_t realSize = inFile.tellg();
    inFile.seekg( 0, std::ios::beg );

    BOOST_CHECK_EQUAL( n * common::typeSize( stype ), realSize );

    inFile.readBinary( data1, n, stype );
    inFile.close();

    BOOST_REQUIRE_EQUAL( data.size(), data1.size() );

    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    AbsType diff = common::Math::real( data.maxDiffNorm( data1 ) );

    BOOST_CHECK_EQUAL( diff, 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
