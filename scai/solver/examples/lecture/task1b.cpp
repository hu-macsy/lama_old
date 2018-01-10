/**
 * @file solver/examples/lecture/task1b.cpp
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
 * @brief ToDo: Missing description in ./solver/examples/lecture/task1b.cpp
 * @author Thomas Brandes
 * @date 15.05.2013
 */

//Solution of task 1b:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

typedef DefaultReal ValueType;

int main()
{
    IndexType size = 4;
    SparseAssemblyStorage<ValueType> sas( size, size, 10 );

    for ( IndexType i = 0; i < size; i++ )
    {
        sas.set( i, i, 2 );
    }

    for ( IndexType i = 0; i < size - 1; i++ )
    {
        sas.set( i + 1, i, 1 );
    }

    for ( IndexType i = 0; i < size - 1; i++ )
    {
        sas.set( i, i + 1, 1 );
    }

    CSRSparseMatrix<ValueType> m ( sas );
    DenseVector<ValueType> rhs( size , 0.0 );
    WriteAccess<ValueType> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; i++ )
    {
        hwarhs[i] = i + 1.0;
    }

    hwarhs.release();
    DenseVector<ValueType> solution( size , 0.0 );
    return 0;
}

