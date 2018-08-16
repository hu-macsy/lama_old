/**
 * @file sparse_dense.cpp
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
 * @brief Benchmark program for performance comparison dense vs sparse vector
 * @author Thomas Brandes
 * @date 09.11.2016
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>

using namespace scai;
using namespace lama;
using namespace std;

using common::Walltime;

/* ----------------------------------------------------------------------------- */

void benchGetRow()
{
    CSRSparseMatrix<double> stencilMatrix;

    MatrixCreator::buildPoisson( stencilMatrix, 3, 27, 30, 30, 30 );

    if ( !stencilMatrix.getColDistribution().isReplicated() )
    {
        dmemo::DistributionPtr repDist( new dmemo::NoDistribution( stencilMatrix.getNumColumns() ) );
        stencilMatrix.redistribute( stencilMatrix.getRowDistributionPtr(), repDist );
    }

    cout << "Poisson matrix: " << stencilMatrix << endl;

    double time1 = Walltime::get();

    DenseVector<double> denseV;

    double res1 = 0;

    { 
        SCAI_REGION( "Main.Dense" )

        for ( IndexType i = 0; i < stencilMatrix.getNumRows(); ++i )
        {
            stencilMatrix.getRow( denseV, i );
            Scalar s = denseV.l2Norm();
            res1 += s.getValue<double>();
        }
    }

    time1 = Walltime::get() - time1;

    double time2 = Walltime::get();

    SparseVector<double> sparseV;

    double res2 = 0;

    { 
        SCAI_REGION( "Main.Sparse" )

        for ( IndexType i = 0; i < stencilMatrix.getNumRows(); ++i )
        {
            stencilMatrix.getRow( sparseV, i );
            Scalar s = sparseV.l2Norm();
            res2 += s.getValue<double>();
        }
    }

    time2 = Walltime::get() - time2;

    cout << "Result = " << res1 << ", " << res2  << endl;

    cout << "Runtime: dense = " << time1 << " seconds, sparse = " << time2 << " seconds" << endl;
}

/* ----------------------------------------------------------------------------- */

void benchGetCol()
{
    ELLSparseMatrix<double> stencilMatrix;

    MatrixCreator::buildPoisson( stencilMatrix, 3, 27, 30, 30, 30 );

    if ( !stencilMatrix.getColDistribution().isReplicated() )
    {
        dmemo::DistributionPtr repDist( new dmemo::NoDistribution( stencilMatrix.getNumColumns() ) );
        stencilMatrix.redistribute( stencilMatrix.getRowDistributionPtr(), repDist );
    }

    cout << "Poisson matrix: " << stencilMatrix << endl;

    double time1 = Walltime::get();

    DenseVector<double> denseV;

    double res1 = 0;

    { 
        SCAI_REGION( "Main.Dense" )

        for ( IndexType i = 0; i < stencilMatrix.getNumRows(); ++i )
        {
            stencilMatrix.getColumn( denseV, i );
            Scalar s = denseV.l2Norm();
            res1 += s.getValue<double>();
        }
    }

    time1 = Walltime::get() - time1;

    double time2 = Walltime::get();

    SparseVector<double> sparseV;

    double res2 = 0;

    { 
        SCAI_REGION( "Main.Sparse" )

        for ( IndexType i = 0; i < stencilMatrix.getNumRows(); ++i )
        {
            stencilMatrix.getColumn( sparseV, i );
            Scalar s = sparseV.l2Norm();
            res2 += s.getValue<double>();
        }
    }

    time2 = Walltime::get() - time2;

    SCAI_ASSERT_EQUAL( res1, res2, "Different results when using dense/sparse vector for getCol" )

    cout << "Runtime: dense = " << time1 << " seconds, sparse = " << time2 << " seconds" << endl;
}

void benchAdd()
{
    const IndexType n = 3000000;

    dmemo::DistributionPtr dist( new dmemo::NoDistribution(  n  ) );

    DenseVector<double> denseV;
    DenseVector<double> denseV1;
    DenseVector<double> denseV2;

    float fillRate = 0.02;

    denseV1.setRandom( dist, fillRate );
    denseV2.setRandom( dist, fillRate );

    SparseVector<double> sparseV;
    SparseVector<double> sparseV1( denseV1 );
    SparseVector<double> sparseV2( denseV2 );

    sparseV = sparseV1 + 2 * sparseV2;
    denseV = denseV1 + 2 * denseV2;

    DenseVector<double> compareV( sparseV );

    SCAI_ASSERT_EQ_ERROR( 0.0, compareV.getLocalValues().maxDiffNorm( denseV.getLocalValues() ), "sparse/dense have different result" );
}

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.bench" )

    common::Settings::parseArgs( argc, argv );

    // evaluate SCAI_NUM_THREADS

    int nThreads;

    if ( common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    if ( false )
    {
        benchGetRow();
    }

    if ( true )
    {
        benchGetCol();
    }

    if ( false )
    {
        benchAdd();
    }
}
