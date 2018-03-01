/**
 * @file solver/examples/lecture/task4.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./solver/examples/lecture/task4.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 4:

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/tracing.hpp>

#include <iostream>

using namespace scai;
using namespace dmemo;
using namespace solver;
using namespace hmemo;
using namespace lama;

typedef DefaultReal ValueType;

int main ( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        return EXIT_FAILURE;
    }

    auto comm = Communicator::getCommunicatorPtr();

    // Read a sparse matrix from the file that has been specified by command line argument

    auto matrix = read<CSRSparseMatrix<ValueType>>( argv[1] );
    std::cout << *comm << ": read matrix: " << matrix << std::endl;

    // Create communicator and distribution
    IndexType size = matrix.getNumRows ( );
    DistributionPtr dist( new BlockDistribution( size, comm ) );
    matrix.redistribute( dist, dist );
    std::cout << *comm << ": redistributed matrix: " << matrix << std::endl;

    // Create solution vector
    auto solution = fill<DenseVector<ValueType>>( dist, 1 );
    std::cout << *comm << ": vector solution : " << solution << std::endl;

    // Compute the rhs that fits our solution to be able to calculate the error later
    auto rhs = eval<DenseVector<ValueType>>( matrix * solution );
    std::cout << *comm << ": vector rhs : " << rhs << std::endl;
    // Forget the solution, i.e. reset solution to zero so that there is something to solve
    solution = ValueType( 0 );

    // Allocate a common logger that prints convergenceHistory
    bool isDisabled = comm->getRank() > 0;
    LoggerPtr logger( new CommonLogger( "CGLogger: ", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly, isDisabled ) );
    // Create a CG solver with logger
    CG<ValueType> cgSolver ( "CGTestSolver", logger );
    // Create a stopping criterion for the iterative solver cgSolver
    NormPtr<ValueType> norm( new L2Norm<ValueType>( ) );

    auto criterion = std::make_shared<ResidualThreshold<ValueType>>( norm, 1E-8, ResidualCheck::Absolute );
    cgSolver.setStoppingCriterion ( criterion );

    // Initialize the solver with the matrix
    cgSolver.initialize ( matrix );
    // Solve, i.e. find solution for given rhs
    cgSolver.solve ( solution, rhs );

    // calculate the error and its L2-Norm
    auto error = fill<DenseVector<ValueType>>( dist, 1 );
    error = error - solution;
    RealType<ValueType> resNorm = l2Norm( error );
    if ( comm->getRank() == 0 )
    {
        std::cout << "L2-Norm of error is " << resNorm << std::endl;
    }

    return EXIT_SUCCESS;
}
