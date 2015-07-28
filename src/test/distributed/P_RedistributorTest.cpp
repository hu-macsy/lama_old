/**
 * @file P_RedistributorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class P_RedistributorTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/TestMacros.hpp>

#include <lama/distribution/Distribution.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/distribution/Redistributor.hpp>

using namespace lama;
using namespace memory;
using namespace common;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

struct P_RedistributorTestConfig
{
    P_RedistributorTestConfig()
    {
        comm = Communicator::get();
    }

    ~P_RedistributorTestConfig()
    {
    }

    CommunicatorPtr comm;
};

BOOST_FIXTURE_TEST_SUITE( P_RedistributorTest, P_RedistributorTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.P_RedistributorTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, ValueType, test_types )
{
    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    IndexType blockLocalSize = distBlock->getLocalSize();
    IndexType cyclicLocalSize = distCyclic->getLocalSize();
    LAMAArray<ValueType> myData1( blockLocalSize );
    {
        WriteAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            data[i] = static_cast<ValueType>( 100 * comm->getRank() + i );
        }
    }
    LAMAArray<ValueType> myData2( cyclicLocalSize );
    Redistributor r1( distCyclic, distBlock );
    Redistributor r2( distBlock, distCyclic );
    LAMA_LOG_DEBUG( logger, "redistribute 1" );
    r1.redistribute( myData2, myData1 );
    LAMA_LOG_DEBUG( logger, "redistribute 2" );
    r2.redistribute( myData1, myData2 );
    {
        ReadAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            ValueType expected = static_cast<ValueType>( 100 * comm->getRank() + i );
            BOOST_CHECK_CLOSE( data[i], expected, 1 );
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
    LAMA_WRITEAT_TEST( r );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
