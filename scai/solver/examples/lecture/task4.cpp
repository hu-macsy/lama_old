/**
 * @file solver/examples/lecture/task4.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task4.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 4:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/tracing.hpp>

#include <iostream>

using namespace scai;
using namespace dmemo;
using namespace solver;
using namespace hmemo;
using namespace lama;

typedef RealType ValueType;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    CSRSparseMatrix<ValueType> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();

    CommunicatorPtr comm( Communicator::getCommunicatorPtr() );
    DistributionPtr dist( new BlockDistribution( size, comm ) );
    m.redistribute( dist, dist );

    DenseVector<ValueType> rhs( size , 0.0 );
    WriteAccess<ValueType> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = ValueType( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    rhs.redistribute( dist );
    DenseVector<ValueType> solution( size, 0.0 );
    solution.redistribute( dist );
    std::cout << "Vector solution : " << solution << std::endl;
    CG cgSolver( "CGTestSolver" );
    CriterionPtr criterion( new IterationCount( 10 ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );

    std::cout << "The solution is: ";

    for ( int i = 0; i < solution.size(); ++i )
    {
        std::cout << solution.getValue( i ) << " ";
    }

    std::cout << std::endl;
    return 0;
}

