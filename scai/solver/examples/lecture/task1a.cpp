/**
 * @file solver/examples/lecture/task1a.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task1a.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 1a:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

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

    //Read a sparse matrix from the passed input file
    CSRSparseMatrix<ValueType> m( argv[1] );

    IndexType size = m.getNumRows();

    DenseVector<ValueType> rhs( size , 0.0 );
    WriteAccess<ValueType> hwarhs( rhs.getLocalValues() );	
    for ( IndexType i = 0; i < size; ++i )
    { 
        hwarhs[i] = i + 1.0;
    }

    hwarhs.release();

    DenseVector<ValueType> solution( size , 0.0 );

    return 0;
}

