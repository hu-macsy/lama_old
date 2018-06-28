/**
 * @file solver.cpp
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
 * @brief Example program for using a solver
 * @author Kai Buschulte
 * @date 12.09.2012
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

// Solver related includes
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/CG.hpp>

#include <iostream>
#include <cstdlib>

using scai::IndexType;
using namespace scai::lama;
using namespace scai::solver;

int main()
{
    const IndexType N1    = 100;
    const IndexType N2    = 100;
    const IndexType NIter = 140;
    //
    // Define the ValueType used for the vector
    //
    typedef double ValueType;
    //
    // Create DenseVectors for solution and right hand side
    //
    DenseVector<ValueType> solution( N1 * N2, ValueType( 2 ) );
    const DenseVector<ValueType> rhs( solution.getDistributionPtr(), ValueType( 1 ) );
    //
    // Create a 2D Poisson Matrix with 9 points and dimension 100 in every direction
    //
    CSRSparseMatrix<ValueType> matrix;
    MatrixCreator::buildPoisson2D( matrix, 9, N1, N2 );
    //
    // Set up CG solver
    //
    // Logger (optional) prints out some solver information to console
    LoggerPtr logger(
        new CommonLogger(
            "<CG>: ",
            LogLevel::solverInformation,
            LoggerWriteBehaviour::toConsoleOnly,
            std::shared_ptr<Timer>( new Timer() ) ) );
    // stopping criterion for solver: stopps after NIter iterations
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( NIter ) );
    // solver itself
    CG<ValueType> solver( "mySolverName", logger );
    // initialization
    solver.setStoppingCriterion( criterion );
    solver.initialize( matrix );
    // solution phase
    solver.solve( solution, rhs );
    solution.writeToFile( "solution.txt" );
    std::cout << "Solution vector is written to 'solution.txt'" << std::endl;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}
