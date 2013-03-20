/**
 * @file InverseSolverTest.cpp
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
 * @brief Contains the implementation of the class InverseSolverTest.
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/InverseSolver.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/expression/VectorExpressions.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<float,double> test_types;

#define LAMA_TO_TEXT( t ) #t

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( InverseSolverTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.InverseSolverTest" );

/* --------------------------------------------------------------------- */

template<typename mt>
void testSolveMethod( ContextPtr context )
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();

    DenseVector<ValueType> solution( system.coefficients.getNumRows(), 0.0 );
    DenseVector<ValueType> reference( system.solution );
    MatrixType coefficients( system.coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "InverseSolverTest uses context = " << context->getType() );

    DenseVector<ValueType> rhs( system.rhs );

    InverseSolver inverseSolver( "InverseSolverTest solver" );
    inverseSolver.initialize( coefficients );

    inverseSolver.solve( solution, rhs );

    DenseVector<ValueType> diff( reference - solution );
    Scalar maxDiff = maxNorm( diff );
    BOOST_CHECK( maxDiff.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( InverseTest, T, test_types ) {
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveMethod< DenseMatrix<ValueType> >( context );
        testSolveMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveMethod< JDSSparseMatrix<ValueType> >( context );
        testSolveMethod< COOSparseMatrix<ValueType> >( context );
        //TODO: DIA crashs
        //testSolveMethod< DIASparseMatrix<ValueType> >( context );

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( InverseTest2, T, test_types ) {
    typedef T ValueType;

    EquationHelper::EquationSystem<ValueType> system =
        EquationHelper::get4x4SystemA<ValueType>();

    const IndexType n = 4;

    DenseVector<ValueType> solution( n, 1.0 );
    DenseVector<ValueType> solution2( n, 1.0 );

    std::string s = "DataType";
    InverseSolver inverseSolver( "InverseSolverTest<" + s + "> solver" );

    DenseMatrix<ValueType> inverse = DenseMatrix<ValueType>( system.coefficients );
    DenseMatrix<ValueType> origin = DenseMatrix<ValueType>( system.coefficients );
    DenseMatrix<ValueType> result = DenseMatrix<ValueType>( system.coefficients );

    inverseSolver.initialize( origin );
    inverseSolver.computeInverse( inverse );

    origin.matrixTimesMatrix( result, 1.0, inverse, 0.0, result );

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            Scalar scalar = result.getValue( i, j );
            if ( i == j )
            {
                BOOST_CHECK_CLOSE( 1.0, scalar.getValue<ValueType>(), 1 );
            }
            else
            {
                BOOST_CHECK( scalar.getValue<ValueType>() < 1E-6 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    InverseSolver inverseSolver( "InverseSolverTest solver" );
    LAMA_WRITEAT_TEST( inverseSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    InverseSolver inverseSolver1( "InverseSolver" );

    SolverPtr solverptr = inverseSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "InverseSolver" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
