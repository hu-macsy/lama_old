/**
 * @file solver/examples/lecture/task0.cpp
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
 * @brief Initial version of a solver to be used as task 0 in lecture.
 * @author Thomas Brandes
 * @date 15.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace lama;
using namespace solver;

typedef DefaultReal ValueType;

int main ( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        return EXIT_FAILURE;
    }

    // Read a sparse matrix from the file that has been specified by command line argument
    auto matrix = read<CSRSparseMatrix<ValueType>>( argv[1] );
    std::cout << "Read matrix: " << matrix << std::endl;

    IndexType size = matrix.getNumRows ( );
    // Create solution vector
    auto solution = fill<DenseVector<ValueType>>( size, 1 );
    std::cout << "Vector solution : " << solution << std::endl;
    // Compute the rhs that fits our solution to be able to calculate the error later
    auto rhs = eval<DenseVector<ValueType>>( matrix * solution );
    std::cout << "Vector rhs : " << rhs << std::endl;
    // Forget the solution, i.e. reset solution to zero so that there is something to solve
    solution = ValueType( 0 );

    // Create a CG solver
    CG<ValueType> cgSolver ( "CGTestSolver" );
    // Create a stopping criterion for the iterative solver cgSolver
    auto norm = std::make_shared<L2Norm<ValueType>>( );
    const ValueType eps = 1E-8;
    auto criterion = std::make_shared<ResidualThreshold<ValueType>>( norm, eps, ResidualCheck::Absolute );
    cgSolver.setStoppingCriterion ( criterion );

    // Initialize the solver with the matrix
    cgSolver.initialize ( matrix );
    // Solve, i.e. find solution for given rhs
    cgSolver.solve ( solution, rhs );

    // calculate the error and its L2-Norm
    auto error = fill<DenseVector<ValueType>>( size, 1 );
    error = error - solution;
    std::cout << "L2-Norm of error is " << l2Norm ( error ) << std::endl;

    return EXIT_SUCCESS;
}

