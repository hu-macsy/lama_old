/**
 * @file LUSolverTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the class LUSolverTest.
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/LUSolver.hpp>
#include <lama/solver/InverseSolver.hpp>

#include <lama/DenseVector.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/Communicator.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/EquationHelper.hpp>

#include <iostream>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LUSolverTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.LUSolverTest" );

/* --------------------------------------------------------------------- */

template<typename T>
void testSolveMethod( EquationHelper::EquationSystem<T> system, const IndexType tilesize )
{
    typedef T ValueType;

    // declare those variable parameters.
    // EquationHelper::EquationSystem<T> system = EquationHelper::get4x4SystemA<DataType>( );
    // const IndexType tilesize = 2;

    const IndexType n = system.coefficients.getNumRows();

    DenseVector<ValueType> solution( n, 1.0 );

    LUSolver luSolver( "LUSolverTest" );
    luSolver.setTileSize( tilesize );

    DenseMatrix<ValueType> coefficients( system.coefficients );

    //BOOST_MESSAGE( "Refactor LUSolver" );
    //luSolver.initialize( coefficients );                                           //TODO: crashes

    DenseVector<ValueType> rhs( system.rhs );

    //  luSolver.solve( solution, rhs );                                               //TODO: crashes

    DenseVector<ValueType> reference( system.solution );

    DenseVector<ValueType> diff( n, 1.0 );

    diff = reference - solution;

    Scalar s = maxNorm( diff );

//    BOOST_CHECK_CLOSE( 0.0, s.getValue<ValueType>(), 1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, T, test_types ) {
    typedef EquationHelper::EquationSystem<T>( *t_system )( );
    const std::vector<t_system>& systems = EquationHelper::getFunctions<T>( );

    for( typename std::vector<t_system>::size_type i=0; i<systems.size( ); ++i )
    {
        const EquationHelper::EquationSystem<T> system = systems[i]( );
        for( int k=1; k<=system.coefficients.getNumRows( ); ++k )
        {
            testSolveMethod( system, 1);
        }
    }
}

/* --------------------------------------------------------------------- */

template<typename T>
void testLUMethod( EquationHelper::EquationSystem<T> system, const IndexType tilesize )
{
    typedef T ValueType;

    // declare those variable parameters.
//    EquationHelper::EquationSystem<T> system = EquationHelper::get4x4SystemA<DataType>( );
//    const IndexType tilesize = 2;

    const IndexType n = system.coefficients.getNumRows();

    DenseVector<ValueType> solution( n, 1.0 );
    DenseVector<ValueType> solution2( n, 1.0 );

    LUSolver luSolver( "LUSolverTest" );
    luSolver.setTileSize( tilesize );

    InverseSolver inverseSolver( "OpenMPInverseSolver" );

    std::vector<IndexType> permutation;

    DenseMatrix<ValueType> result = DenseMatrix<ValueType>( system.coefficients );
    // luSolver.factorMatrixToLU( result, permutation );                                   //TODO: crashes

    DenseMatrix<ValueType> luSolution = DenseMatrix<ValueType>( system.coefficients );
    // ompInverseSolver.computeLUDecomposition( luSolution, permutation );                 //TODO: crashes

    for( IndexType i = 0; i < n; ++i )
    {
        for( IndexType j = 0; j < n; ++j )
        {
            Scalar solScalar = luSolution.getValue( i, j );
            Scalar comScalar = result.getValue( i, j );
            //         BOOST_CHECK_CLOSE( solScalar.getValue<ValueType>(), comScalar.getValue<ValueType>(), 1 );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( luTest, T, test_types) {
    typedef EquationHelper::EquationSystem<T>( *t_system )( );
    const std::vector<t_system>& systems = EquationHelper::getFunctions<T>( );

    for( typename std::vector<t_system>::size_type i=0; i<systems.size( ); ++i )
    {
        const EquationHelper::EquationSystem<T> system = systems[i]( );
        for( int k=1; k<=system.coefficients.getNumRows( ); ++k )
        {
            testLUMethod( system, k );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
