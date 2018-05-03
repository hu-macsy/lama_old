/**
 * @file solver/examples/lecture/task1.cpp
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
 * @brief Solution of lecture task 1.
 * @author Thomas Brandes
 * @date 15.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/VectorAssembly.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace scai::lama;
using namespace scai::solver;

typedef DefaultReal ValueType;

int main( int, char*[] )
{
    const IndexType size = 100;
    const ValueType s    = 0.5;

    MatrixAssembly<ValueType> mAssembly;
        
    for ( IndexType i = 1; i < size - 1; i++ )
    {   
        mAssembly.push( i, i - 1, -s );
        mAssembly.push( i, i, 1 + 2 * s );
        mAssembly.push( i, i + 1, -s );
    }

    mAssembly.push( 0, 0, 1 );
    mAssembly.push( size - 1, size - 1, 1 );

    auto matrix = convert<CSRSparseMatrix<ValueType>>( mAssembly.buildGlobalCOO( size, size ) );

    std::cout << "matrix = " << matrix << std::endl;
    matrix.writeToFile( "matrix.txt" );
    
    const ValueType Tleft  = 1;
    const ValueType Tright = 10;

    VectorAssembly<ValueType> vAssembly;

    vAssembly.push( 0, Tleft );
    vAssembly.push( size - 1, Tright );

    auto rhs = fill<DenseVector<ValueType>>( size , 0 );
    rhs.fillFromAssembly( vAssembly );

    std::cout << "Vector rhs = " << rhs << std::endl;
    rhs.writeToFile( "rhs.txt" );

    auto solution = fill<DenseVector<ValueType>>( size, 0 );

    auto logger = std::make_shared<CommonLogger>( "<CG>: ", LogLevel::solverInformation, LoggerWriteBehaviour::toConsoleOnly );
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

    auto residual = eval<DenseVector<ValueType>>( rhs - matrix * solution );

    std::cout << "Resiudual = " << l2Norm( residual ) << std::endl;

    solution.writeToFile( "solution.txt" );

    std::cout << "Solution has been written to file solution.txt" << std::endl;
}
