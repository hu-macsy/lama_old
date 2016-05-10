/**
 * @file solver/examples/lecture/task2a.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task2a.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 2a:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/tracing.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include <iostream>

using namespace scai::lama;
using namespace scai::solver;
using namespace scai::hmemo;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    CSRSparseMatrix<double> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();
    DenseVector<double> rhs( size , 0.0 );
    WriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    DenseVector<double> solution( size , 0.0 );
    std::cout << "Vector solution : " << solution << std::endl;
    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    CG cgSolver( "CGTestSolver" );
    cgSolver.setStoppingCriterion( rt );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );
    std::cout << "The solution is: ";
    ReadAccess<double> hra( solution.getLocalValues() );

    for ( int i = 0; i < solution.size(); i++ )
    {
        std::cout << hra[i] << " ";
    }

    std::cout << std::endl;

    hra.release();
    return 0;
}

