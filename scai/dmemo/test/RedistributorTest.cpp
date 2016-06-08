/**
 * @file RedistributorTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class RedistributorTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace scai::common;

/* --------------------------------------------------------------------- */

struct RedistributorTestConfig
{
    RedistributorTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
    }

    ~RedistributorTestConfig()
    {
    }

    CommunicatorPtr comm;
};

BOOST_FIXTURE_TEST_SUITE( RedistributorTest, RedistributorTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.RedistributorTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributeTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    IndexType blockLocalSize = distBlock->getLocalSize();
    IndexType cyclicLocalSize = distCyclic->getLocalSize();
    HArray<ValueType> myData1( blockLocalSize );
    {
        WriteAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            data[i] = static_cast<ValueType>( 100 * comm->getRank() + i );
        }
    }
    HArray<ValueType> myData2( cyclicLocalSize );
    Redistributor r1( distCyclic, distBlock );
    Redistributor r2( distBlock, distCyclic );
    SCAI_LOG_DEBUG( logger, "redistribute 1" );
    r1.redistribute( myData2, myData1 );
    SCAI_LOG_DEBUG( logger, "redistribute 2" );
    r2.redistribute( myData1, myData2 );
    {
        ReadAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            ValueType expected = static_cast<ValueType>( 100 * comm->getRank() + i );
            SCAI_CHECK_CLOSE( data[i], expected, 1 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    Redistributor r( distCyclic, distBlock );

    std::ostringstream out;

    out << r ;

    BOOST_CHECK( out.str().length() >  0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
