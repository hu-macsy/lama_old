/**
 * @file lama/examples/labelrank/labelrank.cpp
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
 * @brief ToDo: Missing description in ./lama/examples/labelrank/labelrank.cpp
 * @author Thomas Brandes
 * @date 10.06.2013
 */
#include "Config.hpp"

#include <scai/lama/matrix/all.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/common/Walltime.hpp>

#include <cstdio>
#include <cstdlib>
#include <iomanip>

using namespace std;

using namespace scai;
using namespace lama;
using namespace hmemo;

using common::Walltime;
using common::UnaryOp;

/** ValueType is the type used for matrix and vector elements. */

typedef DefaultReal ValueType;

/** Specify the matrix type for the label matrices.
 *
 *  DenseMatrix is more efficient than CSRSparseMatrix
 */

typedef CSRSparseMatrix<ValueType> LabelMatrix;

/** Specify the matrix type for the affinity matrix.
 *
 *  CSRSparseMatrix is more efficient than DenseMatrix.
 */

typedef CSRSparseMatrix<ValueType> AffinityMatrix;

/** Parameters of this inference algorithm */

struct Properties
{
    size_t logLevel;  //!<  control amount of output
    double threshold; //!< termination criterion
} props;

/** Determine evidence indexes of Dense label matrix. */
/* UNUSED
static void getEvidence( std::vector<IndexType>& evidenceVector, const DenseStorage<ValueType>& labelsMatrix )
{
    // Get indexes of all rows where at least one element is found

    evidenceVector.clear();

    const IndexType numRows = labelsMatrix.getNumRows();
    const IndexType numCols = labelsMatrix.getNumColumns();

    ReadAccess<ValueType> dense( labelsMatrix.getData() );

    // if labelsMatrix[i,:] has entry, we add row index i

    for ( IndexType i = 0; i < numRows; ++i )
    {
        bool evidence = false;

        for ( int j = 0; j < numCols; ++j )
        {
            if ( dense[ i * numCols + j ] == 1 )
            {
                // cout << "evidence at row " << i;
                // cout << " j = " << j ;
                // cout << ", val = " << dense[ i * numCols + j] << endl;
                evidence = true;
                break;
            }
        }

        if ( evidence )
        {
            evidenceVector.push_back( i );
        }
    }

    cout << "evidence vector, len = " << evidenceVector.size() << endl;
}
*/

/** Determine evidence indexes of CSR label matrix. */

static void getEvidence( std::vector<IndexType>& evidenceVector, const CSRStorage<ValueType>& labelsMatrix )
{
    // Get indexes of all rows where at least one element is found
    evidenceVector.clear();
    const IndexType numRows = labelsMatrix.getNumRows();
    //const IndexType numCols = labelsMatrix.getNumColumns();
    ReadAccess<IndexType> csrIA( labelsMatrix.getIA() );

    // if labelsMatrix[i,:] has entry, we add row index i

    for ( IndexType i = 0; i < numRows; ++i )
    {
        bool evidence = csrIA[i] < csrIA[i + 1];

        if ( evidence )
        {
            evidenceVector.push_back( i );
        }
    }

    cout << "evidence vector, len = " << evidenceVector.size() << endl;
}

/** Update affinity matrix if it is dense. */
/* UNUSED
static void update( DenseStorage<ValueType>& affinityMatrix, const std::vector<IndexType>& evidenceVector )
{
    // This is some kind of hack to make sure that labelsMatrix is not updated in
    // all rows where at least one element is found

    //const IndexType numRows = affinityMatrix.getNumRows();
    const IndexType numCols = affinityMatrix.getNumColumns();

    WriteAccess<ValueType> affinityDense( affinityMatrix.getData() );

    // if labelsMatrix[i,:] has entry, we set affinityMatrix[i,:] = [0, ..., 1 ... ]
    // so labelsMatrix[i,:] remains unchanged for affinityMatrix * labelsMatrix

    IndexType size = evidenceVector.size();

    for ( IndexType ii = 0; ii < size; ++ii )
    {
        IndexType i = evidenceVector[ii];

        // set row i of affinityMatrix to identity

        for ( IndexType j = 0; j < numCols; ++j )
        {
            affinityDense[ i * numCols + j ] = 0.0;
        }

        affinityDense[ i * numCols + i ] = 1.0;
    }
}
*/

/** Update affinity matrix if it is CSR. */

static void update( CSRStorage<ValueType>& affinityMatrix, const std::vector<IndexType>& evidenceVector )
{
    // This is some kind of hack to make sure that labelsMatrix is not updated in
    // all rows where at least one element is found

    const IndexType size = evidenceVector.size();

    hmemo::HArray<IndexType> ia;
    hmemo::HArray<IndexType> ja;
    hmemo::HArray<ValueType> values;

    IndexType numRows;
    IndexType numColumns;

    affinityMatrix.splitUp( numRows, numColumns, ia, ja, values );

    {
        ReadAccess<IndexType> csrIA( ia );
        WriteAccess<IndexType> csrJA( ja );
        WriteAccess<ValueType> csrValues( values );

        // if labelsMatrix[i,:] has entry, we set affinityMatrix[i,:] = [0, ..., 1 ... ]
        // so labelsMatrix[i,:] remains unchanged for affinityMatrix * labelsMatrix

        for ( IndexType ii = 0; ii < size; ++ii )
        {
            IndexType i = evidenceVector[ii];

            // set row i of affinityMatrix to identity

            // SCAI_ASSERT_ERROR( csrIA[i] < csrIA[i+1], "row " << i << " of affinity matrix empty" );

            if ( csrIA[i] >= csrIA[i + 1] )
            {
                continue;
            }

            csrJA[ csrIA[i] ] = i;
            csrValues[ csrIA[i] ] = 1.0;

            for ( IndexType jj = csrIA[i] + 1; jj < csrIA[i + 1]; ++jj )
            {
                csrJA[ jj ] = i;
                csrValues[ jj ] = 0.0;
            }
        }
    }

    affinityMatrix = CSRStorage<ValueType>( numRows, numColumns, std::move( ia ), std::move( ja ), std::move( values ) );
}

int main( int argc, char* argv[] )
{
    // config
    Config config;
    // set properties
    props.logLevel = 20;
    props.threshold = 1e-5;
    const char* wFilename;
    const char* yFilename;

    if ( argc >= 3 )
    {
        wFilename = argv[1];
        yFilename = argv[2];

        for ( int i = 3; i < argc; ++i )
        {
            config.setArg( argv[i] );
        }
    }
    else
    {
        cout << "Usage: " << argv[0] << " <affinity.mtx> <labels.mtx>" << endl;
        exit( 1 );
    }

    double start = Walltime::get();

    auto affinityMatrix = read<AffinityMatrix>( wFilename );

    cout << "loading affinityMatrix took " << Walltime::get() - start << " secs." << endl;
    cout << "affinityMatrix = " << affinityMatrix << endl;

    start = Walltime::get();

    auto labelsMatrix = read<LabelMatrix>( yFilename );

    cout << "loading labelsMatrix took " << Walltime::get() - start << " secs." << endl;
    cout << "labelsMatrix = " << labelsMatrix << endl;
    // compute diagonal degree matrix and invert it
    start = Walltime::get();
    //const IndexType numRows = affinityMatrix.getNumRows();
    const IndexType numCols = affinityMatrix.getNumColumns();

    auto oneVector = denseVectorFill<ValueType>( numCols, 1.0 );
    auto y         = denseVectorEval( affinityMatrix * oneVector );  // rowSums

    y = 1 / y; 

    affinityMatrix.scaleRows( y );  // scales each row
    cout << "invert/scale calculations took " << Walltime::get() - start << " secs." << endl;
    // update affinityMatrix so that labelsMatrix remains unchanged
    // where initial entries are set
    std::vector<IndexType> evidences;
    getEvidence( evidences, labelsMatrix.getLocalStorage() );
    cout << "#evidences = " << evidences.size() << endl;
    update( affinityMatrix.getLocalStorage(), evidences );
    cout << "Starting Label Propagation..." << endl;
    LabelMatrix labelsMatrixNew;

    if ( props.logLevel < 30 )
    {
        cout << "| " << "iter " << " | " << "max_diff" << " | " << "nnz      " << " | " << "nnz      " << " | " << "secs.   " << " | " << endl;
        cout << "-------------------------------------------------------" << endl;
    }

    int iteration = 0;
    cout << "Running on " << config.getContext() << endl;
    // do all operations on matrices at CUDA device where possible
    labelsMatrix.setContextPtr( config.getContextPtr() );
    labelsMatrixNew.setContextPtr( config.getContextPtr() );
    affinityMatrix.setContextPtr( config.getContextPtr() );
    double totalStart = Walltime::get();

    while ( true )
    {
        iteration += 1;
        double start = Walltime::get();
        // (D−1 * W) * Y(t)
        // w already stores (D-1 * W)
        labelsMatrixNew = affinityMatrix * labelsMatrix;
        int nnzOld = labelsMatrix.getNumValues();
        ValueType maxDiff = labelsMatrix.maxDiffNorm( labelsMatrixNew );
        // use more efficient swap instead of: labelsMatrix = labelsMatrixNew;
        // Note: works only for same matrix type
        labelsMatrix.swap( labelsMatrixNew );

        if ( props.logLevel < 30 )
        {
            cout << "| ";
            cout << std::setw( 5 ) << iteration << " | ";
            cout << std::setiosflags( std::ios::fixed ) << std::setprecision( 6 ) << std::setw( 8 ) << maxDiff << " | ";
            cout << std::setw( 9 ) << nnzOld << " | ";
            int nnz = labelsMatrix.getNumValues();
            cout << std::setw( 9 ) << nnz << " | ";
            cout << std::setiosflags( std::ios::fixed ) << std::setprecision( 6 ) << std::setw( 8 ) << ( Walltime::get() - start ) << " | ";
            cout << endl;
        }

        if ( maxDiff < props.threshold )
        {
            if ( props.logLevel < 30 )
            {
                cout << "-------------------------------------------------------" << endl;
            }

            cout << "Converged after " << iteration << " iterations with max diff " << maxDiff << endl;
            break;
        }

        if ( iteration == config.mMaxIters )
        {
            if ( props.logLevel < 30 )
            {
                cout << "-------------------------------------------------------" << endl;
            }

            cout << "Stopped without convergence" << endl;
            break;
        }
    }

    cout << "saving new labels to hdd..." << endl;
    labelsMatrix.writeToFile( "labels-new.mtx" );
    cout << "total run time was " << Walltime::get() - totalStart << " secs." << endl;
}
