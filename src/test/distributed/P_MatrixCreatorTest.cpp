/**
 * @file P_MatrixCreatorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class P_MatrixCreatorTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 06.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/TestMacros.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/matutils/MatrixCreator.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<double,float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_MatrixCreatorTestConfig
{
    P_MatrixCreatorTestConfig()
    {
        comm = CommunicatorFactory::get();
    }

    ~P_MatrixCreatorTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_MatrixCreatorTest, P_MatrixCreatorTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.P_MatrixCreatorTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, T, test_types ) {
    typedef T ValueType;

    const IndexType globalSize = 100;

    DistributionPtr dist( new BlockDistribution( globalSize, comm ) );

    CSRSparseMatrix<ValueType> matrix;

    matrix.allocate( dist, dist );

    MatrixCreator<ValueType>::fillRandom( matrix, 0.05 );

//  std::cout << "random matrix = " << matrix << std::endl;
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testPoisson, T, test_types ) {
    typedef T ValueType;

    CSRSparseMatrix<ValueType> matrix;

    LAMA_CHECK_THROW(
    {   MatrixCreator<T>::buildPoisson1D( matrix, 5, 100 );}, Exception );

    MatrixCreator<ValueType>::buildPoisson1D( matrix, 3, 100 );

//std::cout << "Poisson1D3P matrix = " << matrix << std::endl;

    MatrixCreator<ValueType>::buildPoisson2D( matrix, 5, 10, 10 );

//std::cout << "Poisson2D5P matrix = " << matrix << std::endl;

    MatrixCreator<ValueType>::buildPoisson3D( matrix, 7, 5, 5, 5 );

//std::cout << "Poisson3D7P matrix = " << matrix << std::endl;
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
