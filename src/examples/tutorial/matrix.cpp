/**
 * @file matrix.cpp
 *
 * @license
 * Copyright (c) 2012
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief matrix.cpp is an example matrix vector multiplication (CSR matrix format).
 * @author Bea Hornef
 * @date 16.05.2013
 * $Id$
 */

#include <lama.hpp>

// Matrix & vector related includes
#include <lama/DenseVector.hpp>
#include <lama/expression/all.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/storage/CSRStorage.hpp>

using namespace lama;

int main()
{
	// EXAMPLE multiplication of a dense vector with a sparse matrix in CSR format.

	typedef float ValueType;
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
    ValueType vectorValues[] =
    {   6.0f, 4.0f, 7.0f, -9.3f};

    // All data has to be stored in LAMAArrays.
    const LAMAArray<IndexType> matrixIA = LAMAArray<IndexType>( numRows + 1, ia );
    const LAMAArray<IndexType> matrixJA = LAMAArray<IndexType>( numValues, ja );
    const LAMAArray<ValueType> mValues  = LAMAArray<ValueType>( numValues, matrixValues );
    const LAMAArray<ValueType> vValues  = LAMAArray<ValueType>( numColumns, vectorValues );

    //  Alternative code for the next 2 lines:
    //  CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( numRows, numColumns, numValues, matrixIA, matrixJA, matrixValues );

    // Create a CSRStorage.
    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>();
    csrStorage->setCSRData( numRows, numColumns, numValues, matrixIA, matrixJA, mValues );

    // Allocate and fill vector for the multiplication.
    DenseVector<ValueType> vector( numColumns, 0.0 );
    vector.setValues(vValues);
    // Allocation of the result vector.
    DenseVector<ValueType> result( numRows, 0.0);

    // Distribution pointer are needed to construct a CSRSparseMatrix.
    lama::DistributionPtr rowDist(new lama::NoDistribution(numRows));
    lama::DistributionPtr colDist(new lama::NoDistribution(numColumns));

    // Allocation of the CSRSparseMatrix.
    CSRSparseMatrix<ValueType> csrMatrix(*csrStorage, rowDist, colDist);

    // The multiplication itself.
    result = csrMatrix * vector;
    result.writeToFile( "result.txt" , File::FORMATTED );
}
