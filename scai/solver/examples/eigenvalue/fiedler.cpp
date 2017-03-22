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

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/CG.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// import common 
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <stdlib.h>

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
void makeLaplacian( CSRSparseMatrix<double>& L )
{
    using namespace hmemo;

    // set all non-zero non-diagonal elements to -1
    // set all diagonal elements to number of non-diagonal elements

    CSRStorage<double>& localStorage = L.getLocalStorage();

    const HArray<IndexType>& ia = localStorage.getIA();
    const HArray<IndexType>& ja = localStorage.getJA();
    HArray<double>& values = localStorage.getValues();

    {
        const double minusOne = -1;

        ReadAccess<IndexType> rIA( ia );
        ReadAccess<IndexType> rJA( ja );
        WriteAccess<double> wValues( values );

        for ( IndexType i = 0; i < localStorage.getNumRows() ; ++i )
        {
            const double sum = rIA[i+1] - rIA[i] - 1;
 
            for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
            {
                wValues[jj] = rJA[jj] == i ? sum : minusOne;
            }
        }
    }

    // check that x = (1, 1, ..., 1 ) is eigenvector with eigenvalue 0

    DenseVector<double> x( L.getColDistributionPtr(), 1 );
    DenseVector<double> y( L * x );

    SCAI_ASSERT_LT_ERROR( y.maxNorm(), Scalar( 1e-8 ), "L not Laplacian matrix" )
}

/** CG solver where L * w is replaced with L * w - s' * w * r - r' * w * s */

void solve( Vector& x, Matrix& L, Vector& r, Vector& s, const Vector& b, const Scalar& eps )
{
    DenseVector<double> res ( b - L * x );

    res += s.dotProduct( x ) * r;
    res += r.dotProduct( x ) * s;
    res[0] = 0.0;

    DenseVector<double> d ( res );

    Scalar rOld = res.dotProduct( res );
    L2Norm norm; 
    
    IndexType maxIter = 20;

    Scalar rNorm = norm( res );

    for ( int k = 0 ; k < maxIter and rNorm > eps; k++ )
    {
        DenseVector<double> z( L * d );

        z -= s.dotProduct( d ) * r;
        z -= r.dotProduct( d ) * s;
        z[0] = 0.0;

        Scalar alpha = rOld / d.dotProduct( z );
        x = x + alpha * d;
        res = res - alpha * z;
        Scalar rNew = res.dotProduct( res );
        Scalar beta = rNew / rOld;
        d = res + beta * d;
        rOld = rNew; 
        rNorm = norm( res );
        std::cout << "Iter k = " << k << " : norm( r ) = " << rNorm.getValue<double>() << std::endl;
    }
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
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    SCAI_ASSERT_GE_ERROR( argc, 2, "insufficient arguments, call " << argv[0] << " <filename>" )

    std::string filename = argv[1];

    CSRSparseMatrix<double> L;  // laplacian matrix

    IndexType kmax = 100;        // maximal number of iterations
    Scalar    eps  = 1e-7;       // accuracy for maxNorm 

    L.readFromFile( argv[1] );

    SCAI_ASSERT_EQ_ERROR( L.getNumRows(), L.getNumColumns(), "matrix not square" )

    makeLaplacian( L );   // make the matrix Laplacian and check it

    const IndexType n = L.getNumRows();

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr blockDist( new dmemo::BlockDistribution( n, comm ) );
    L.redistribute( blockDist, blockDist );

    DenseVector<double> u( L.getRowDistributionPtr(), 1 );

    double n12 = common::Math::sqrt( double( n  ) );

    u[0] = n12 + 1;
 
    Scalar alpha = n + n12;

    DenseVector<double> h ( L * u );
    h /= alpha;

    Scalar gamma = u.dotProduct( h ) / alpha * 0.5;

    DenseVector<double> v( h - gamma * u );

    DenseVector<double> r( u ); r[0] = 0.0;
    DenseVector<double> s( v ); s[0] = 0.0;

    DenseVector<double> t( L.getRowDistributionPtr(), 1.0 );

    DenseVector<double> y;
    DenseVector<double> diff;

    t[0] = 0.0;
    DenseVector<double> z( t );
    
    for ( IndexType k = 0; k < kmax; ++k )
    {
        // normalize t

        t = t / t.l2Norm();

        y = L * t;
        y[0] = 0.0;                  // fill element as we actually use L[2:n,2:n] 
        y -= s.dotProduct( t ) * r;  
        y -= r.dotProduct( t ) * s;

        Scalar lambda = t.dotProduct( y );
        diff = y - lambda * t;
        Scalar diffNorm = diff.maxNorm();

        std::cout << "Iter " << k << ", lambda = " << lambda << ", diff = " << diffNorm << std::endl;

        if ( diffNorm < eps )
        {
            break;
        }

        solve( z, L, r, s, t, eps );

        t = z;
    }

    t[0] = 0.0;
    Scalar beta = u.dotProduct( t ) / alpha;
    t = t - beta * u;

    t.writeToFile( "eigenvector.mtx" );
}
