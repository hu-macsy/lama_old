/**
 * @file GlobalAddressingPlanTest.cpp
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
 * @brief Test routines for the class GlobalAddressingPlan
 * @author Thomas Brandes
 * @date 12.12.2018
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/data/test_case.hpp>

#include <scai/dmemo/GlobalAddressingPlan.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/test/TestMacros.hpp>


using namespace scai;
using namespace dmemo;
using namespace hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GlobalAddressingPlanTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GlobalAddressingPlanTest" )

/* --------------------------------------------------------------------- */

/** 
 *   Help function to create a distributed array
 *
 *   @param[in] dist is the distribution of the array
 *   @param[in] fill is a function that returns for a global index the value
 *   @returns   the local part owned by this processor of the distributed array
 */
template<typename ValueType>
static HArray<ValueType> distributedArray( const Distribution& dist, ValueType ( *fill )( IndexType ) )
{
    HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {
        IndexType localIndex = 0;   // running local index

        for ( auto& entry : hostWriteOnlyAccess( localArray, dist.getLocalSize() ) )
        {
            entry = fill( dist.local2Global( localIndex++ ) );
        }

    }  // filled the local array with 'global' values

    return localArray;    // each processor gets its local part
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( gatherTest )
{
    typedef DefaultReal ValueType;

    const IndexType N  = 100;
    const IndexType m  = 15;

    auto dist = blockDistribution( N );

    HArray<IndexType> indexes( m );
    utilskernel::HArrayUtils::setRandom( indexes, N - 1 );

    HArray<ValueType> localArray;

    auto fillArray = []( IndexType k ) { return ValueType( 2 * k + 1 ); };

    HArray<ValueType> remoteArray = distributedArray<ValueType>( *dist, fillArray );

    HArray<ValueType> expArray;  // will contain expected gathered values

    {
        auto expected = hostWriteOnlyAccess( expArray, m );

        auto rIndexes = hostReadAccess( indexes );

        IndexType i = 0;

        for ( const auto index : hostReadAccess( indexes )  )
        {
            expected[ i++ ] = fillArray( index );
        }
    }

    auto plan = globalAddressingPlan( *dist, indexes );

    localArray.resize( m );
    plan.gather( localArray, remoteArray );

    BOOST_TEST( hostReadAccess( localArray ) == hostReadAccess( expArray ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scatterTest )
{
    typedef DefaultReal ValueType;

    const IndexType N  = 100;

    auto dist = blockDistribution( N );
  
    // use another distribution to have a global permutation of indexes

    auto dist1 = cyclicDistribution( N, 2 );
    auto indexes = dist1->ownedGlobalIndexes();

    // define a fill function for the distributed array

    auto fillArray = []( IndexType k ) { return ValueType( 2 * k + 1 ); };

    HArray<ValueType> localArray = distributedArray<ValueType>( *dist1, fillArray );
    HArray<ValueType> expArray   = distributedArray<ValueType>( *dist, fillArray );

    auto plan = globalAddressingPlan( *dist, indexes );

    HArray<ValueType> remoteArray( dist->getLocalSize(), ValueType( 0 ) );

    plan.scatter( remoteArray, localArray, common::BinaryOp::COPY );

    BOOST_TEST( hostReadAccess( remoteArray ) == hostReadAccess( expArray ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
