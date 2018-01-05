/**
 * @file solver/examples/lecture/task2.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task2.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 2:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/all.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/tracing.hpp>

#include <iostream>

using namespace scai::lama;
using namespace scai::hmemo;

typedef RealType ValueType;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    int maxIter = 10;   // maximal number of iterations

    if ( argc > 2 )
    {
        sscanf( argv[2], "%d", &maxIter );
    }

    CSRSparseMatrix<ValueType> A( argv[1] );
    std::cout << "Read matrix A : " << A << std::endl;
    IndexType size = A.getNumRows();
    DenseVector<ValueType> b( size, 0 );
    {
        WriteAccess<ValueType> writeB( b.getLocalValues() );

        for ( IndexType i = 0; i < size; ++i )
        {
            // writeB[i] = 1;
            writeB[i] = ValueType( i + 1 );
        }
    }
    std::cout << "Vector b : " << b << std::endl;
    DenseVector<ValueType> x( size , 0.0 );
    std::cout << "Vector x : " << x << std::endl;
    // d = r = b - A * x
    // help = A * x;
    DenseVector<ValueType> r ( b - A * x );
    DenseVector<ValueType> d ( r );
    ValueType rOld = r.dotProduct( r );
    ValueType eps = 0.00001;
    L2Norm<ValueType> norm;

    for ( int k = 0 ; k < maxIter and norm( r ) > eps; k++ )
    {
        DenseVector<ValueType> z( A * d );
        ValueType alpha = rOld / d.dotProduct( z );
        x = x + alpha * d;
        r = r - alpha * z;
        ValueType rNew = r.dotProduct( r );
        ValueType beta = rNew / rOld;
        d = r + beta * d;
        rOld = rNew;
        ValueType rnorm = norm( r );
        std::cout << "Iter k = " << k << " : norm( r ) = " << rnorm << std::endl;
    }

    return 0;
}

