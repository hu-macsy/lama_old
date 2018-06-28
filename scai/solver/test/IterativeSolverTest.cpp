/**
 * @file IterativeSolverTest.cpp
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
 * @brief Generic tests for solvers derived from IterativeSolver.
 * @author Jan Ecker
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/TypeTraits.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/IterativeSolver.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <memory>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace lama;
using namespace solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( IterativeSolverTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.IterativeSolverTest" )

// ---------------------------------------------------------------------------------------------------------------

/** This routine provides a 2-dim stencil matrix for all tests */

template<typename ValueType>
CSRSparseMatrix<ValueType> setupMatrix( const IndexType N1, const IndexType N2 )
{
    CSRSparseMatrix<ValueType> matrix;                    // gets default context
    MatrixCreator::buildPoisson2D( matrix, 9, N1, N2 );
    return matrix;
}

// ---------------------------------------------------------------------------------------------------------------

/** Help class where an object is a vector of all available iterative solvers. */

template<typename ValueType>
class IterativeSolvers : public std::vector<SolverPtr<ValueType> >
{

public:

    /** Constructor creates vector of all available solvers. */

    IterativeSolvers()
    {
        std::vector<SolverCreateKeyType> values;

        _Solver::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); i++ )
        {
            if ( values[i].second == "Kaczmarz" )
            {
                continue;   
            }

            _SolverPtr _solver( _Solver::create( values[i] ) );

            std::shared_ptr<IterativeSolver<ValueType> > itSolver =
               std::dynamic_pointer_cast<IterativeSolver<ValueType> >( _solver );

            if ( itSolver )
            {
                this->push_back( itSolver );
            }
        }
    }

    // Destructor will free all matrix storages due to use of shared pointers

    IterativeSolver<ValueType>& get( const size_t i )
    {
        SCAI_ASSERT_LT( i, this->size(), "Index out of range" )
        SolverPtr<ValueType> solver = this->at( i );
        IterativeSolver<ValueType>* iterativeSolver = dynamic_cast<IterativeSolver<ValueType>*>( solver.get() );
        SCAI_ASSERT( iterativeSolver, "Not iterative solver" )
        return *iterativeSolver;
    }
};

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( DefaultCriterionTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const CSRSparseMatrix<ValueType> matrix = setupMatrix<ValueType>( 40, 40 );

    SCAI_LOG_INFO( logger, "DefaultCriterionTest with matrix = " << matrix )

    // setup must deliver a square matrix with same distribution for row and columns

    BOOST_REQUIRE_EQUAL( matrix.getRowDistribution(), matrix.getColDistribution() );

    auto rhs = fill<DenseVector<ValueType>>( matrix.getColDistributionPtr(), 1, matrix.getContextPtr() );
    auto solution = fill<DenseVector<ValueType>>( matrix.getRowDistributionPtr(), 0, matrix.getContextPtr() );

    solution.setContextPtr( matrix.getContextPtr() );

    IterativeSolvers<ValueType> solvers;  // All available iterative solvers

    for ( size_t i = 0; i < solvers.size(); i++ )
    {
        IterativeSolver<ValueType>& iterativeSolver = solvers.get( i );

        SCAI_LOG_INFO( logger, "Testing solver " << iterativeSolver );

        iterativeSolver.initialize( matrix );
        iterativeSolver.solve( solution, rhs );

        BOOST_CHECK_EQUAL( iterativeSolver.getIterationCount(), IndexType( 1 ) );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( IterationCountStoppingCriterionTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const CSRSparseMatrix<ValueType> matrix = setupMatrix<ValueType>( 40, 40 );

    SCAI_LOG_INFO( logger, "matrix = " << matrix );

    const ValueType solutionInitValue = 1.0;

    auto solution      = fill<DenseVector<ValueType>>( matrix.getColDistributionPtr(), solutionInitValue, matrix.getContextPtr() );
    auto exactSolution = eval<DenseVector<ValueType>>( solution + 1 );
    auto rhs           = eval<DenseVector<ValueType>>( matrix * exactSolution );

    IndexType numIterations = 300;
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( numIterations ) );

    IterativeSolvers<ValueType> solvers;  // vector with incarnation of each solver

    for ( size_t i = 0; i < solvers.size(); i++ )
    {
        IterativeSolver<ValueType>& iterativeSolver = solvers.get( i );

        solution = solutionInitValue;   // reset solution for each solver

        SCAI_LOG_INFO( logger, "Testing solver " << iterativeSolver );

        // TODO: this should be tested for ALL solvers (not only iterative)
        // Solver has to be initialized before solve is called
        //BOOST_CHECK_THROW ( {iterativeSolver->solve( solution, rhs );}, scai::common::Exception );

        iterativeSolver.setStoppingCriterion( criterion );
        iterativeSolver.initialize( matrix );
        iterativeSolver.solve( solution, rhs );

        BOOST_CHECK_EQUAL( numIterations, iterativeSolver.getIterationCount() );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConservativeTest )
{
    // Idea:  solver.solve( exactSolution, rhs ), 1 step, exactSolution must be the same

    typedef SCAI_TEST_TYPE ValueType;

    const CSRSparseMatrix<ValueType> matrix = setupMatrix<ValueType>( 40, 40 );

    SCAI_LOG_INFO( logger, "matrix = " << matrix );

    ValueType solutionValue = 1;

    auto exactSolution = fill<DenseVector<ValueType>>( matrix.getColDistributionPtr(), solutionValue, matrix.getContextPtr() );
    auto solution      = DenseVector<ValueType>( exactSolution );
    auto rhs           = eval<DenseVector<ValueType>>( matrix * exactSolution );

    IndexType numIterations = 1;
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( numIterations ) );

    IterativeSolvers<ValueType> solvers;  // vector with incarnation of each solver

    for ( size_t i = 0; i < solvers.size(); i++ )
    {
        IterativeSolver<ValueType>& iterativeSolver = solvers.get( i );

        SCAI_LOG_INFO( logger, "ConservativeTest: solver " << iterativeSolver );

        iterativeSolver.setStoppingCriterion( criterion );

        solution = solutionValue + scai::common::TypeTraits<ValueType>::eps1();

        iterativeSolver.initialize( matrix );
        iterativeSolver.solve( solution, rhs );

        ValueType maxDiff = maxNorm( iterativeSolver.getResidual() );

        SCAI_LOG_INFO( logger, "ConservativeTest: solver " << iterativeSolver << ", maxDiff = " << maxDiff )

        ValueType eps = 1e-4;

        BOOST_CHECK_MESSAGE( maxDiff < eps,
                             "Solver: " << iterativeSolver << " max diff = " << maxDiff << " >= " << eps );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( SolveTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const CSRSparseMatrix<ValueType> matrix = setupMatrix<ValueType>( 10, 10 );

    SCAI_LOG_INFO( logger, "matrix = " << matrix );

    const ValueType solutionInitValue = 1;

    auto colDist       = matrix.getColDistributionPtr();
    auto solution      = fill<DenseVector<ValueType>>( colDist, solutionInitValue, matrix.getContextPtr() );
    auto exactSolution = fill<DenseVector<ValueType>>( colDist, solutionInitValue + 1, matrix.getContextPtr() );
    auto rhs           = eval<DenseVector<ValueType>>( matrix * exactSolution );

    IndexType maxExpectedIterations = 3000;

    auto criterion = std::make_shared<IterationCount<ValueType>>( maxExpectedIterations );

    IterativeSolvers<ValueType> solvers;  // vector with incarnation of each solver

    for ( size_t i = 0; i < solvers.size(); i++ )
    {
        IterativeSolver<ValueType>& iterativeSolver = solvers.get( i );

        SCAI_LOG_INFO( logger, "Testing solver " << iterativeSolver );

        iterativeSolver.setStoppingCriterion( criterion );
        iterativeSolver.initialize( matrix );
        solution = solutionInitValue;
        iterativeSolver.solve( solution, rhs );
        auto diff = eval<DenseVector<ValueType>>( solution - exactSolution );
        RealType<ValueType> realMaxNorm     = maxNorm( diff );
        RealType<ValueType> expectedMaxNorm = 1E-4;
        SCAI_LOG_INFO( logger, "maxNorm of diff = ( solution - exactSolution ) = " << realMaxNorm );
        BOOST_CHECK_MESSAGE( realMaxNorm < expectedMaxNorm,
                             "Solver: " << iterativeSolver << " max norm = " << realMaxNorm << " >= " << expectedMaxNorm );
        // Test solve with preconditioner
        SolverPtr<ValueType> preconditioner( new TrivialPreconditioner<ValueType>( "Trivial preconditioner" ) );
        iterativeSolver.setPreconditioner( preconditioner );
        // TODO: this should be tested for ALL preconditioners
        // Solver has to be initialized AFTER a preconditioner is set
        //BOOST_CHECK_THROW ( {iterativeSolver.solve( solution, rhs );}, scai::common::Exception );
        iterativeSolver.initialize( matrix );
        solution = solutionInitValue;
        iterativeSolver.solve( solution, rhs );
        BOOST_CHECK_EQUAL( maxExpectedIterations, iterativeSolver.getIterationCount() );
        diff = solution - exactSolution;
        realMaxNorm     = maxNorm( diff );
        expectedMaxNorm = 1E-4;
        SCAI_LOG_INFO( logger, "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << realMaxNorm );
 
        BOOST_CHECK_MESSAGE( realMaxNorm < expectedMaxNorm,
                             "Solver: " << iterativeSolver << " max norm = " << realMaxNorm << " >= " << expectedMaxNorm );
    }
}

BOOST_AUTO_TEST_SUITE_END();
