/**
 * @file matrix.cpp
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
 * @brief matrix.cpp is an example matrix vector multiplication (CSR matrix format).
 * @author Bea Hornef
 * @date 16.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace hmemo;
using namespace lama;

//
// EXAMPLE multiplication of a dense vector with a sparse matrix in CSR format.
//
int main()
{
    //
    // Define the ValueType used for the vector
    //
    typedef DefaultReal ValueType;

    // initialize matrix and vector values
    //
    //    matrix( 7 x 4 )   vector( 4 )
    //     6  0   0  4          6
    //     7  0   0  0          4
    //     0  0  -9  4          7
    //     2  5   0  3         -9
    //     2  0   0  1 
    //     0  0   0  0 
    //     0  1   0  2

    IndexType numRows    =  7;
    IndexType numColumns =  4;
    IndexType numValues  = 12;

    // VERY IMPORTANT for CSR format: ia is an offset array adding number of values for each row starting with 0.
    // ja stores the column position of each value within a row

    IndexType rawIA[]     = {   0,    2,    3,    5,       8,   10, 10,    12};
    IndexType rawJA[]     = {   0, 3, 0,    2, 3, 0, 1, 3, 0, 3,    1,  3 };
    ValueType rawValues[] = {   6, 4, 7,   -9, 4, 2, 5, 3, 2, 1,    1,  2 };

    // Vector values for our multiplication.
    ValueType rawVector[] = {   6, 4, 7, -9 };

    // All raw data is copied into heterogeneous arrays to support operations on any device

    HArray<IndexType> csrIA( numRows + 1, rawIA );
    HArray<IndexType> csrJA( numValues, rawJA );
    HArray<ValueType> csrValues( numValues, rawValues );
    HArray<ValueType> vectorValues( numColumns, rawVector );

    // Create a CSRStorage.

    CSRStorage<ValueType> csrStorage ( numRows, numColumns, csrIA, csrJA, csrValues );

    // Create the CSRSparseMatrix.

    CSRSparseMatrix<ValueType> csrMatrix( csrStorage );

    // Allocate and fill vector for the multiplication.

    DenseVector<ValueType> vector( vectorValues );

    // Allocation of the result vector.

    DenseVector<ValueType> result( numRows, 0 );

    //
    // The multiplication itself.
    //
    result = csrMatrix * vector;
    result = 2 * transpose( csrMatrix ) * vector + 2 * result;
    //
    // print vector to file result.txt (simple ASCII format)
    //
    result.writeToFile( "result.txt" );

    std::cout << "DenseVector result values have been written to 'result.txt'" << std::endl;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}
