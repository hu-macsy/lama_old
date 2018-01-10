/**
 * @file solver/examples/lecture/task3.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task3.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 3:

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>

#include <iostream>

using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace solver;

typedef DefaultReal ValueType;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    //Read a sparse matrix from the passed input file
    CSRSparseMatrix<ValueType> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();
    //Create rhs vector
    DenseVector<ValueType> rhs( size , 0.0 );
    std::cout << "Vector rhs : " << rhs << std::endl;
    //Create solution vector
    DenseVector<ValueType> solution( size , 1.0 );
    std::cout << "Vector solution : " << solution << std::endl;
    //Compute the rhs that fits our solution to be able to calculate the error later
    rhs = m * solution;
    //Reset solution to zero so that there is something to solve
    solution = 0.0;
    std::shared_ptr<Timer> timer( new Timer() );
    LoggerPtr logger( new CommonLogger( "CGLogger", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly, timer ) );
    //Create a CG solver using the Constructor with logger as an argument
    CG<ValueType> cgSolver( "CGTestSolver", logger );
    //Create a stopping criterion for the iterative solver cgSolver
    NormPtr<ValueType> norm( new L2Norm<ValueType>() );
    CriterionPtr<ValueType> criterion( new ResidualThreshold<ValueType>( norm, 1E-8, ResidualCheck::Absolute ) );
    cgSolver.setStoppingCriterion( criterion );
    //Initialize the solver
    cgSolver.initialize( m );
    //Solve m * solution = rhs
    cgSolver.solve( solution, rhs );
    //calculate the error and its L2-Norm
    DenseVector<ValueType> error( size, 1.0 );
    error = error - solution;
    std::cout << "L2-Norm of error is " << l2Norm( error ) << std::endl;
    return 0;
}

