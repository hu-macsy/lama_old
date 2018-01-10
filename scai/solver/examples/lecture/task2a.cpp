/**
 * @file solver/examples/lecture/task2a.cpp
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

using namespace scai;
using namespace lama;
using namespace solver;
using namespace hmemo;

typedef DefaultReal ValueType;

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
    DenseVector<ValueType> rhs( size , 0.0 );
    WriteAccess<ValueType> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = ValueType( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    DenseVector<ValueType> solution( size , 0.0 );
    std::cout << "Vector solution : " << solution << std::endl;
    ValueType eps = 0.00001;
    NormPtr<ValueType> norm( new L2Norm<ValueType>() );
    CriterionPtr<ValueType> rt( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Absolute ) );
    CG<ValueType> cgSolver( "CGTestSolver" );
    cgSolver.setStoppingCriterion( rt );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );
    std::cout << "The solution is: ";
    ReadAccess<ValueType> hra( solution.getLocalValues() );

    for ( IndexType i = 0; i < solution.size(); i++ )
    {
        std::cout << hra[i] << " ";
    }

    std::cout << std::endl;
    hra.release();
    return 0;
}

