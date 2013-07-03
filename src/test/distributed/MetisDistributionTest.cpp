/**
 * @file MetisDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief MetisDistributionTest.cpp
 * @author Lauretta Schubert
 * @date 01.07.2013
 * since 1.1.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/MetisDistribution.hpp>

#include <test/TestMacros.hpp>
#include <test/Configuration.hpp>
#include <test/distributed/DistributionTest.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <vector>

using namespace lama;
using namespace boost;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<CSRSparseMatrix<float>,ELLSparseMatrix<float>,COOSparseMatrix<float>,JDSSparseMatrix<float>,
        DIASparseMatrix<float>,CSRSparseMatrix<double>,ELLSparseMatrix<double>,COOSparseMatrix<double>,
        JDSSparseMatrix<double>,DIASparseMatrix<double> > MatrixTypes;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct MetisDistributionTestConfig
{
    MetisDistributionTestConfig()
    {
        comm = CommunicatorFactory::get( "MPI" );

        rank = comm->getRank();
        size = comm->getSize();

        std::string prefix = Configuration::getInstance().getPath();
        std::string formattedInputFile = prefix + "/bcspwr01.mtx";
        matrix = CSRSparseMatrix<double>( formattedInputFile );

        globalSize = matrix.getNumRows();

        // weights
        float weight = 1.0 / size;
        parts.reserve( size );
        for ( int i = 0; i < size - 1; ++i )
        {
            parts[i] = weight;
        }
        parts[ size - 1 ] = 1.0 - (size - 1) * weight;

        dist = DistributionPtr( new MetisDistribution<double>( comm, matrix, parts ) );
    }

    ~MetisDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType globalSize;

    DistributionPtr dist;

    CSRSparseMatrix<double> matrix;

    std::vector<float> parts;
};

BOOST_FIXTURE_TEST_SUITE( MetisDistributionTest, MetisDistributionTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.MetisDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run test method " << testcase << " in MetisDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( isEqualTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::ValueType ValueType;

    MatrixType distMatrix( matrix );

    DistributionPtr generaldist1( new MetisDistribution<ValueType>( comm, distMatrix, parts ) );
    DistributionPtr generaldist2( generaldist1 );
    DistributionPtr generaldist3( new MetisDistribution<ValueType>( comm, distMatrix, parts ) );

    BOOST_CHECK(  (*generaldist1).isEqual( *generaldist2 ) );
    BOOST_CHECK( !(*generaldist1).isEqual( *generaldist3 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
