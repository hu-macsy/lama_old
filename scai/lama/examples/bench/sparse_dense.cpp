/**
 * @file sparse_dense.cpp
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

    CSRSparseMatrix<double> stencilMatrix;

    MatrixCreator::buildPoisson( stencilMatrix, 3, 27, 30, 30, 30 );

    // make sure that columns are replicated, otherwise getRow is inefficient at all

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
            stencilMatrix.getRow1( sparseV, i );
            Scalar s = sparseV.l2Norm();
            res2 += s.getValue<double>();
        }
    }

    time2 = Walltime::get() - time2;

    cout << "Result = " << res1 << ", " << res2  << endl;

    cout << "Runtime: dense = " << time1 << " seconds, sparse = " << time2 << " seconds" << endl;
}
