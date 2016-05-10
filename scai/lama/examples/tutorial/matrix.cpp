/**
 * @file matrix.cpp
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
 * @brief matrix.cpp is an example matrix vector multiplication (CSR matrix format).
 * @author Bea Hornef
 * @date 16.05.2013
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::hmemo;
using namespace scai::lama;
using scai::utilskernel::LArray;

//
// EXAMPLE multiplication of a dense vector with a sparse matrix in CSR format.
//
int main()
{
    //
    // Define the ValueType used for the vector
    //
	typedef float ValueType;

	//
	// initialize matrix and vector values
	//

	//    Our matrix:
	//     6, 0,     0, 4,
	//     7, 0,     0, 0,
	//     0, 0, -9.3f, 4,
	//     2, 5,     0, 3,
	//     2, 0,     0, 1,
	//     0, 0,     0, 0,
	//     0, 1,     0, 2

    IndexType numRows    = 7;
	IndexType numColumns = 4;
    IndexType numValues  = 12;

    // VERY IMPORTANT for CSR format: ia is an offset array adding values of each row starting with 0.
    // ja stores the position of each value within a line,
    IndexType ia[] =
    {   0, 2, 3, 5, 8, 10, 10, 12};
    IndexType ja[] =
    {   0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3};
    ValueType matrixValues[] =
    {   6.0f, 4.0f, 7.0f, -9.3f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };

    // Vector values for our multiplication.
    ValueType vectorValues[] = {   6.0f, 4.0f, 7.0f, -9.3f };

    // All data has to be stored in LAMA Arrays.
    const LArray<IndexType> matrixIA = LArray<IndexType>( numRows + 1, ia );
    const LArray<IndexType> matrixJA = LArray<IndexType>( numValues, ja );
    const LArray<ValueType> mValues  = LArray<ValueType>( numValues, matrixValues );
    const LArray<ValueType> vValues  = LArray<ValueType>( numColumns, vectorValues );

    // Create a CSRStorage.
    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( numRows, numColumns, numValues,
                                                                   matrixIA, matrixJA, mValues );

    //  Alternative code for the last line
//    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>();
//    csrStorage->setCSRData( numRows, numColumns, numValues, matrixIA, matrixJA, mValues );

    // Allocate and fill vector for the multiplication.
    DenseVector<ValueType> vector( numColumns, 0.0 );
    vector.setValues( vValues );
    // Allocation of the result vector.
    DenseVector<ValueType> result( numRows, 0.0 );

    // Distribution pointer are needed to construct a CSRSparseMatrix.
    scai::dmemo::DistributionPtr rowDist( new scai::dmemo::NoDistribution( numRows ) );
    scai::dmemo::DistributionPtr colDist( new scai::dmemo::NoDistribution( numColumns ) );

    // Allocation of the CSRSparseMatrix.
    CSRSparseMatrix<ValueType> csrMatrix( *csrStorage, rowDist, colDist );

    //
    // The multiplication itself.
    //
    result = csrMatrix * vector;

    //
    // print vector to file result.frm/.vec (SAMG format)
    //
    result.writeToFile( "result" , File::FORMATTED );

    std::cout << "DenseVector is written to 'result.frm/.vec'" << std::endl;

    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

    return EXIT_SUCCESS;
}
