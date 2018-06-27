/**
 * @file fiedler.cpp
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
 * @brief Inverse power method incorporated with Householder deflation
 * @author Thomas Brandes
 * @date 22.03.2017
 */

#include <scai/lama.hpp>

#include <scai/solver/examples/eigenvalue/HouseholderTransformedMatrix.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/CG.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// we use own instrumentation with SCAI_REGION here
#include <scai/tracing.hpp>

// import common 
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <stdlib.h>

typedef scai::DefaultReal ValueType;  // is double, only float if double not instantiated

using namespace scai;
using namespace lama;
using namespace solver;

/** This routine converts any symmetric matrix with diagonal elements
 *  to a laplacian matrix.
 *
 *  \code
 *     0.5  0.1  -  -        1   -1   -   -
 *     0.2  1.3  1.4 -       -1   2   1   - 
 *      -   0.3  1.1 1.2      -  -1   2   1
 *      -    -   -0.1 0.2    -   -   -1   2 
 *  \endcode
 */
void makeLaplacian( CSRSparseMatrix<ValueType>& L )
{
    SCAI_REGION( "main.makeLaplacian" )

    using namespace hmemo;
    using namespace utilskernel;

    // set all non-zero non-diagonal elements to -1
    // set all diagonal elements to number of non-diagonal elements

    CSRStorage<ValueType>& localStorage = L.getLocalStorage();

    const HArray<IndexType>& ia = localStorage.getIA();
    const HArray<IndexType>& ja = localStorage.getJA();

    // Okay, not very safe but we know what we do, and certainly do not change the size 

    HArray<ValueType>& values = const_cast<HArray<ValueType>&>( localStorage.getValues() );

    {
        const ValueType minusOne = -1;

        ReadAccess<IndexType> rIA( ia );
        ReadAccess<IndexType> rJA( ja );
        WriteAccess<ValueType> wValues( values );

        for ( IndexType i = 0; i < localStorage.getNumRows() ; ++i )
        {
            const ValueType sum = static_cast<ValueType>( rIA[i+1] - rIA[i] - 1 );
 
            for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
            {
                wValues[jj] = rJA[jj] == i ? sum : minusOne;
            }
        }
    }

    // check that x = (1, 1, ..., 1 ) is eigenvector with eigenvalue 0

    auto x = fill<DenseVector<ValueType>>( L.getColDistributionPtr(), 1 );
    auto y = eval<DenseVector<ValueType>>( L * x );

    SCAI_ASSERT_LT_ERROR( y.maxNorm(), ValueType( 1e-8 ), "L not Laplacian matrix" )
}

/** Main program to determine the Fiedler vector for a Laplacian matrix
 *
 *  Method: Inverse power method incorporated with Householder deflation
 *
 *  Paper: An efficient and accurate method to compute the Fiedler vector based on 
 *         Householder deflation and inverse power iteration
 *         Jian-ping Wu, Jun-qiang Song, Wei-min Zhang 
 */
int main( int argc, const char* argv[] )
{
    SCAI_REGION( "main.Fiedler" )

    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    SCAI_ASSERT_GE_ERROR( argc, 2, "insufficient arguments, call " << argv[0] << " <filename>" )

    std::string filename = argv[1];

    CSRSparseMatrix<ValueType> L;  // laplacian matrix

    IndexType kmax = 100;        // maximal number of iterations

    auto eps = common::Math::real( ValueType( 1e-5 ) );  // accuracy for maxNorm, must not be complex

    L.readFromFile( argv[1] );

    double time = common::Walltime::get();

    SCAI_ASSERT_EQ_ERROR( L.getNumRows(), L.getNumColumns(), "matrix not square" )

    makeLaplacian( L );   // make the matrix Laplacian and check it

    const auto n = L.getNumRows();

    auto blockDist = std::make_shared<dmemo::BlockDistribution>( n );

    L.redistribute( blockDist, blockDist );

    auto u = fill<DenseVector<ValueType>>( L.getRowDistributionPtr(), 1 ); 

    ValueType n12 = common::Math::sqrt( ValueType( n  ) );

    u[0] = n12 + 1;
 
    ValueType alpha = n + n12;

    HouseholderTransformedMatrix<ValueType> HLH( L, u, alpha );

    auto t = fill<DenseVector<ValueType>>( L.getRowDistributionPtr(), 1.0 );

    DenseVector<ValueType> y;
    DenseVector<ValueType> diff;

    t[0] = 0.0;

    DenseVector<ValueType> z( t );
    
    // set up the CG solver that is used for the inverse power method

    CG<ValueType> cgSolver( "InversePowerMethodSolver" );

    auto criterion1 = std::make_shared<IterationCount<ValueType>>( 50 );
    NormPtr<ValueType> norm( Norm<ValueType>::create( "L2" ) );           // Norm from factory
    auto criterion2 = std::make_shared<ResidualThreshold<ValueType>>( norm, eps, ResidualCheck::Absolute );
    auto criterion  = std::make_shared<Criterion<ValueType>>( criterion1, criterion2, BooleanOp::OR );

    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( HLH );

    for ( IndexType k = 0; k < kmax; ++k )
    {
        SCAI_REGION( "main.iterate" )

        // normalize t

        t = t / t.l2Norm();

        y = HLH * t;

        auto lambda = t.dotProduct( y );

        diff = y - lambda * t;

        auto diffNorm = diff.maxNorm();

        std::cout << "Iter " << k << ", lambda = " << lambda << ", diff = " << diffNorm << std::endl;

        if ( diffNorm < eps )
        {
            break;
        }

        // solve( z, HLH, t, eps );
        cgSolver.solve( z, t );

        t = z;
    }

    t[0] = 0.0;
    ValueType beta = u.dotProduct( t ) / alpha;
    t = t - beta * u;

    time = common::Walltime::get() - time ;

    std::cout << "Time: " << time << " seconds" << std::endl;

    t.writeToFile( "eigenvector.mtx" );
}
