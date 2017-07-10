/**
 * @file sparsekernel/test/BuildSparseTest.cpp
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
    BuildSparseVector<ValueType> myVector( 10 );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );

    myVector.push( 1, 5 );

    BOOST_CHECK_EQUAL( IndexType( 1 ), myVector.getLength() );
    BOOST_CHECK( !myVector.isEmpty() );

    IndexType i;
    ValueType v;

    myVector.pop( i, v );

    BOOST_CHECK_EQUAL( static_cast<IndexType>( 1 ), i );
    BOOST_CHECK_EQUAL( static_cast<ValueType>( 5 ), v );

    BOOST_CHECK_EQUAL( IndexType( 0 ), myVector.getLength() );
    BOOST_CHECK( myVector.isEmpty() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
