/**
 * @file lstest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Test program for least square solver 
 * @author Thomas Brandes
 * @date 05.03.2017
 * @since 1.0.0
 */

#include <scai/solver/examples/leastsquare/GramianMatrix.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/CGNE.hpp>
#include <scai/solver/CGNR.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include <iostream>

using namespace scai;
using namespace lama;
using namespace solver;
using namespace dmemo;

int main( int, char** ) 

{
    typedef double ValueType;

    // coefficients of the matrix A

    ValueType rawA[] = { 0.9501,   0.7620,   0.6153,   0.4057,
                         0.2311,   0.4564,   0.7919,   0.9354,
                         0.6068,   0.0185,   0.9218,   0.9169,
                         0.4859,   0.8214,   0.7382,   0.4102,
                         0.8912,   0.4447,   0.1762,   0.8936 };

    // coefficients of the rhs b

    ValueType rawB[] = { 0.0578, 0.3528, 0.8131, 0.0098, 0.1388 };

    // expected solution

    ValueType rawX[] = { 0.1887, -0.7166, 0.5604, 0.2107 };

    IndexType numRows = 5;
    IndexType numColumns = 4;

    CSRSparseMatrix<ValueType> m;

    hmemo::HArrayRef<double> data( numRows * numColumns, rawA );

    DistributionPtr rowDist( new NoDistribution( numRows ) );
    DistributionPtr colDist( new NoDistribution( numColumns ) );

    m.setDenseData( rowDist, colDist, data, 0.0 );

    GramianMatrix mTm( m );

    DenseVector<ValueType> b( numRows, rawB );
    DenseVector<ValueType> b1( b * m );

    DenseVector<ValueType> x0( colDist, 0 );
    DenseVector<ValueType> bestX( numColumns, rawX );


    // definie stopping criteria

    double eps = 1e-8;
    CriterionPtr criterion1( new IterationCount( 20 ) );
    NormPtr norm( Norm::create( "L2" ) );   // Norm from factory
    CriterionPtr criterion2( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    CriterionPtr criterion( new Criterion( criterion1, criterion2, Criterion::OR ) );

    // define common logger to print convergence history

    LoggerPtr logger( new CommonLogger ( "",
                                         LogLevel::convergenceHistory,
                                         LoggerWriteBehaviour::toConsoleOnly
                                       ) );

    // Do it with CG

    CG solver( "NormalEquationSolver" );
    solver.setLogger( logger );
    solver.setStoppingCriterion( criterion );
    solver.initialize( mTm );

    solver.solve( x0, b1 );

    DenseVector<ValueType> computedRes( m * x0 - b );
    DenseVector<ValueType> expectedRes( m * bestX - b );

    std::cout << "residual norm of computed x : " << computedRes.l2Norm() << std::endl;
    std::cout << "residual norm of expected x : " << expectedRes.l2Norm() << std::endl;

    // Do it with CGNR

    CGNR solver1( "NormalEquationSolver" );

    solver1.setStoppingCriterion( criterion );
    solver1.setLogger( logger );

    x0 = 0;  
    solver1.initialize( m );
    solver1.solve( x0, b );

    computedRes =  m * x0 - b;

    std::cout << "residual norm of computed x : " << computedRes.l2Norm() << std::endl;

    return 0;
}
