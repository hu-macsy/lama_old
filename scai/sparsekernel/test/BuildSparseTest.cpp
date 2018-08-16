/**
 * @file sparsekernel/test/BuildSparseTest.cpp
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
 * @brief Contains tests for the classes BuildSparseIndexes and BuildSparseVector
 * @author Thomas Brandes
 * @date 22.12.2016
 */

// boost

#include <boost/test/unit_test.hpp>

// others

#include <scai/sparsekernel/openmp/BuildSparseIndexes.hpp>

#include <scai/logging.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace sparsekernel;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BuildSparseTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BuildSparseTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildSparseIndexesTest )
{
    BuildSparseIndexes myIndexes( 10 );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myIndexes.getLength() );
    BOOST_CHECK( myIndexes.isEmpty() );

    myIndexes.pushIndex( 3 );

    BOOST_CHECK_EQUAL( IndexType( 1 ), myIndexes.getLength() );
    BOOST_CHECK( !myIndexes.isEmpty() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( buildSparseVectorTest, ValueType, scai_numeric_test_types )
{
    BuildSparseVector<ValueType> myVector( 10, ValueType( 1 ) );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );

    myVector.push( 1, 5, common::BinaryOp::MULT );
    myVector.push( 1, 3, common::BinaryOp::MULT );

    BOOST_CHECK_EQUAL( IndexType( 1 ), myVector.getLength() );
    BOOST_CHECK( !myVector.isEmpty() );

    IndexType i;
    ValueType v;

    myVector.pop( i, v );

    BOOST_CHECK_EQUAL( IndexType( 1 ), i );
    BOOST_CHECK_EQUAL( ValueType( 15 ), v );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( buildSparseVectorTest1, ValueType, scai_numeric_test_types )
{
    BuildSparseVector<ValueType> myVector( 10 );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );

    myVector.push( 1, 5, common::BinaryOp::COPY );
    myVector.push( 1, 3, common::BinaryOp::COPY );

    BOOST_CHECK_EQUAL( IndexType( 1 ), myVector.getLength() );
    BOOST_CHECK( !myVector.isEmpty() );

    IndexType i;
    ValueType v;

    myVector.pop( i, v );

    BOOST_CHECK_EQUAL( static_cast<IndexType>( 1 ), i );
    BOOST_CHECK_EQUAL( static_cast<ValueType>( 3 ), v );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
