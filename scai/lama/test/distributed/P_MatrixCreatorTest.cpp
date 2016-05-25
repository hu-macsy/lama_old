/**
 * @file P_MatrixCreatorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Contains the implementation of the class P_MatrixCreatorTest.
 * @author Alexander Büchel, Lauretta Schubert
 * @date 06.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using scai::common::Exception;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_MatrixCreatorTestConfig
{
    P_MatrixCreatorTestConfig()
    {
        comm = Communicator::getCommunicator();
    }

    ~P_MatrixCreatorTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_MatrixCreatorTest, P_MatrixCreatorTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_MatrixCreatorTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, test_types )
{
    const IndexType globalSize = 100;
    DistributionPtr dist( new BlockDistribution( globalSize, comm ) );
    CSRSparseMatrix<ValueType> matrix;
    matrix.allocate( dist, dist );
    MatrixCreator<ValueType>::fillRandom( matrix, 0.05 );
//  std::cout << "random matrix = " << matrix << std::endl;
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testPoisson, ValueType, test_types )
{
    CSRSparseMatrix<ValueType> matrix;
    SCAI_CHECK_THROW(
    {   MatrixCreator<ValueType>::buildPoisson1D( matrix, 5, 100 );}, Exception );
    MatrixCreator<ValueType>::buildPoisson1D( matrix, 3, 100 );
//std::cout << "Poisson1D3P matrix = " << matrix << std::endl;
    MatrixCreator<ValueType>::buildPoisson2D( matrix, 5, 10, 10 );
//std::cout << "Poisson2D5P matrix = " << matrix << std::endl;
    MatrixCreator<ValueType>::buildPoisson3D( matrix, 7, 5, 5, 5 );
//std::cout << "Poisson3D7P matrix = " << matrix << std::endl;
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
