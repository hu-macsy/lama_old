/**
 * @file solver/examples/lecture/task2.cpp
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
 * @brief Implementation of a CG solver using LAMA textbook syntax
 * @author Thomas Brandes
 * @date 15.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/tracing.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace scai::lama;
using namespace scai::hmemo;

typedef DefaultReal ValueType;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        return EXIT_FAILURE;
    }

    int maxIter = 10;   // maximal number of iterations

    if ( argc > 2 )
    {
        maxIter = std::stoi( argv[2], nullptr, 0 );
    }

    auto matrix = read<CSRSparseMatrix<ValueType>>( argv[1] );
    std::cout << "Read matrix : " << matrix << std::endl;
    IndexType size = matrix.getNumRows();
    auto rhs = linearDenseVector<ValueType>( size, 1, 1 );
    std::cout << "Vector rhs : " << rhs << std::endl;
    auto solution = fill<DenseVector<ValueType>>( size, 0 );

    auto r = eval<DenseVector<ValueType>>( rhs - matrix * solution );
    DenseVector<ValueType> d( r );

    ValueType rOld = r.dotProduct( r );

    L2Norm<ValueType> norm;
    RealType<ValueType> eps = 0.00001;
    auto rnorm = norm( r );

    DenseVector<ValueType> z;  // temporary vector

    for ( int k = 0 ; k < maxIter and rnorm > eps; k++ )
    {
        z = matrix * d;
        ValueType alpha = rOld / d.dotProduct( z );
        solution += alpha * d;
        r -= alpha * z;
        ValueType rNew = r.dotProduct( r );
        ValueType beta = rNew / rOld;
        d = r + beta * d;
        rOld = rNew;
        rnorm = norm( r );  // alternatively: common::Math::sqrt( rNew )
        std::cout << "Iter k = " << k << " : norm( r ) = " << rnorm << ", " << std::endl;
    }

    solution.writeToFile( "solution.txt" );

    return EXIT_SUCCESS;
}

