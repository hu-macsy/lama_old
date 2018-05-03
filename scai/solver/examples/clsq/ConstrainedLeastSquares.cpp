/**
 * @file ConstrainedLeastSquares.cpp
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
 * @brief Implementation of methods for class to solve LeastSquares with constraints
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#include <scai/lama/matrix/AbstractMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/common/Math.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/MINRES.hpp>
#include <scai/solver/GMRES.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include <scai/tracing.hpp>

#include <iostream>

#include "ConstrainedLeastSquares.hpp"

namespace scai
{

namespace solver
{ 

using namespace lama;

SCAI_LOG_DEF_LOGGER( ConstrainedLeastSquares::logger, "ConstrainedLeastSquares" )

// return: lb < x < ub, elementwise

bool isInterior( const DenseVector<double>& x,
                 const DenseVector<double>& lb,
                 const DenseVector<double>& ub )
{
    return x.all( common::CompareOp::LT, ub ) && x.all( common::CompareOp::GT, lb );
}

void ConstrainedLeastSquares::dualityGap( 

    double& gap, 
    double& dualObj, 
    const DenseVector<double>& b, 
    const DenseVector<double>& x, 
    const DenseVector<double>& lb, 
    const DenseVector<double>& ub )

{
    SCAI_REGION( "lsqb.dualityGap" )

    // Given x*, compute an estimate for the duality gap by implictly and
    // analytically computing a dual feasible point which converges to the 
    // optimum as x* converges to the optimum for the primal problem

    DenseVector<double> x_s( x - lb );
    DenseVector<double> b_s( b - mA * lb );
    DenseVector<double> u_s( ub - lb  );
    DenseVector<double> res( mA * x_s - b_s );
    DenseVector<double> kappa( 2 * res );

    DenseVector<double> mu;

    if ( mAT.get() )
    {
         mu = ( -1 ) * (*mAT) * kappa;
    }
    else
    {
         mu = ( -1 ) * kappa * mA;
    }

    mu.setScalar( Scalar( 0 ), common::BinaryOp::MAX );

    Scalar sDualObj = - kappa.dotProduct( kappa ) * 0.25  - kappa.dotProduct( b_s ) - mu.dotProduct( u_s );
    Scalar sGap = res.dotProduct( res ) - sDualObj;

    dualObj = sDualObj.getValue<double>();
    gap = sGap.getValue<double>();
}

double ConstrainedLeastSquares::centralPathObjective (
    const DenseVector<double>& b,
    const DenseVector<double>& x,
    const double tau,
    const DenseVector<double>& lb,
    const DenseVector<double>& ub )
{
    SCAI_REGION( "lsqb.centralPathObjective" )

    DenseVector<double> tmp ( x - lb );
    tmp.log(); 
    double s1 = tmp.sum().getValue<double>();
    tmp = ub - x;
    tmp.log(); 
    double s2 = tmp.sum().getValue<double>();
    double barrier = -s1 - s2;
    DenseVector<double> res( mA * x - b );
    double dp = res.dotProduct( res ).getValue<double>();
    double value = tau * dp + barrier;

    std::cout << "central path, t = " << tau << ", resnorm = " << common::Math::sqrt( dp ) 
              << ", barrier = " << barrier << ", value = " << value << std::endl;

    return value;
}

double ConstrainedLeastSquares::stepSize( 
    const DenseVector<double>& b,
    const DenseVector<double>& x,
    const double tau,
    const DenseVector<double>& lb,
    const DenseVector<double>& ub,
    const DenseVector<double>& dx,
    const double alpha,
    const double beta )
{
    SCAI_REGION( "lsqb.stepSize" )

//      function [ s ] = step_size (A, b, Vector x, t, l, u, Vector dx, alpha, beta)

// % Gradient g = g(x, t)
// d = 1 ./ (u - x) - 1 ./ (x - l);
// g = 2 * t * A' * (A * x - b) + d;

    DenseVector<double> d1( ub - x ); d1.invert();
    DenseVector<double> d2( x - lb ); d2.invert();
    DenseVector<double> d( d1 - d2 );

    DenseVector<double> g( mA * x - b );

    if ( mAT.get() )
    {
        g = 2 * tau * (*mAT) * g  + d;
    }
    else
    {
        g = 2 * tau * g * mA + d;
    }

    IndexType k = 0;

    double objValueAtX = centralPathObjective( b, x, tau, lb, ub );

    Scalar objDiff = alpha * g.dotProduct( dx );

    std::cout << "objValue@X = " << objValueAtX << ", diff = " << objDiff << std::endl;

    double s = 1.0;

    while ( true )
    {
        DenseVector<double> x_new( x + s * dx );

        if ( isInterior( x_new, lb, ub ) )
        {
            double objValueAtXNew = centralPathObjective(  b, x_new, tau, lb, ub );

            double objDiff = objValueAtXNew - objValueAtX;

            std::cout << "objValue@Xnew = " << objValueAtXNew << ", diff = " << objDiff << std::endl;

            if ( objValueAtXNew < objValueAtX + alpha * s * g.dotProduct( dx) ) 
            {
                return s;
            }
        }
        else
        {
            std::cout << "x_new not interior for s = " << s << std::endl;
        }

        s *= beta;   // make s smaller

        k = k + 1;

        std::cout << "stepSize( k = " << k << " ) = " << s << std::endl;
    }
}

void ConstrainedLeastSquares::computeSearchDirection(
     DenseVector<double>& dx,
     const DenseVector<double>& b,
     const DenseVector<double>& x,
     const double tau,
     const DenseVector<double>& lb,
     const DenseVector<double>& ub,
     const double gap )
{
    SCAI_REGION ( "clsq.computeSearchDirection" )

    // Gradient g = g(x, t)

    // d = 1 ./ (u - x) - 1 ./ (x - l);

    DenseVector<double> d1( ub - x );  
    d1.invert();
    DenseVector<double> d2( x - lb ); 
    d2.invert();
    DenseVector<double> d( d1 - d2 );

    // g = 2 * t * A' * (A * x - b) + d;

    DenseVector<double> tmp( mA * x - b );

    DenseVector<double> g;

    if ( mAT.get() )
    {
        // use matrix times vector, is faster

        g =  2 * tau * ( *mAT ) * tmp  + d;
    }
    else
    {
        // no explicit transposed matrix, use vector times matrix

        g =  2 * tau * tmp * mA + d;
    }

    g *= -1;

    // Hessian H = H(x, t)
    // D = 1 ./ ((x - l).^2) + 1 ./ ((u - x).^2);
    // H = @(x) 2 * t * A' * (A * x) + D .* x;

    d1 *= d1;
    d2 *= d2;
    DenseVector<double> D( d1 + d2 );

    mCentralPathHessian.update( D, tau );

   // % Preconditioner for H (note that it is diagonal)
   // P_inv = 1 ./ (2 * t * diagATA + D);

   // % Solve H dx = -g
   // % e is just a multiplier for the tolerance
   // e = 0.01;
   // tol = max(1e-14, min(0.1, e * gap / norm(g)));

   /// relative tolerance, 

   // [dx, ~, ~, ~] = pcg(H, - g, tol, 1000, @(x) P_inv .* x);

    // definie stopping criteria

    double eps = 1e-14;
    double e = 0.01;

    double normG = g.l2Norm().getValue<double>();
   
    double tol = common::Math::max( eps, common::Math::min( 0.1, e * gap / normG ) );

    SCAI_LOG_DEBUG( logger, "Use relative tol = " << tol << " for CG" )

    CriterionPtr criterion1( new IterationCount( 1000 ) );
    NormPtr norm( Norm::create( "L2" ) );   // Norm from factory
    CriterionPtr criterion2( new ResidualThreshold( norm, tol, ResidualThreshold::Relative ) );
    CriterionPtr criterion( new Criterion( criterion1, criterion2, Criterion::OR ) );

    // define common logger to print convergence history

    DenseVector<double> diagonal( 2 * tau * mDiagATA + D );

    mDiagonalMatrix.setDiagonal( diagonal );

    std::shared_ptr<Jacobi> preconditioner( new Jacobi( "JacobiPreconditioner" ) );

    preconditioner->initialize( mDiagonalMatrix );

    // Do it with CG

    CG solver( "searchDirectionSolver" );
    solver.setLogger( mSolverLogger );
    solver.setStoppingCriterion( criterion );
    solver.setPreconditioner( preconditioner );
    solver.initialize( mCentralPathHessian );
    solver.solve( dx, g );
}

ConstrainedLeastSquares::ConstrainedLeastSquares( const lama::CSRSparseMatrix<double>& A ) :

    mA( A ),
    mCentralPathHessian( A ),
    mTolerance( 0.01 ), 
    mMaxIter( 1000 )
{
    SCAI_LOG_INFO( logger, "Constructed solver for constrained least squares, A = " << A )

    mSolverLogger.reset( new CommonLogger ( "", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly) );

    // allocate, initialize the matrix used in the preconditioner

    mDiagonalMatrix.setIdentity( A.getColDistributionPtr() );
    mDiagonalMatrix.setCommunicationKind( A.getCommunicationKind() );
    mDiagonalMatrix.setContextPtr( A.getContextPtr() );
}

void ConstrainedLeastSquares::setTolerance( double tolerance )
{
    SCAI_ASSERT_GT_ERROR( tolerance, double( 0 ), " must be positive" )
    mTolerance = tolerance;
}

void ConstrainedLeastSquares::setMaxIter( IndexType maxIter )
{
    mMaxIter = maxIter;
}

void ConstrainedLeastSquares::useTranspose()
{
    SCAI_REGION ( "clsq.useTranspose" )

    SCAI_ASSERT_ERROR( mAT.get() == NULL, "transposed matrix already available" )

    mAT.reset( new CSRSparseMatrix<double>() );
    mAT->setCommunicationKind( mA.getCommunicationKind() );
    mAT->setContextPtr( mA.getContextPtr() );
    mAT->assignTranspose( mA );

    // the central path Hessian matrix might also used the transposed matrix

    mCentralPathHessian.setTransposed( *mAT );
}

void ConstrainedLeastSquares::solve(
    DenseVector<double>& x,
    const DenseVector<double>& b,
    const DenseVector<double>& lb,
    const DenseVector<double>& ub )
{
    SCAI_REGION ( "clsq.solve" )

    const IndexType n = mA.getNumColumns();

    SCAI_ASSERT_EQ_ERROR( n, lb.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, ub.size(), "size mismatch" )

    SCAI_ASSERT_EQ_ERROR( b.size(), mA.getNumRows(), "size mismatch" )

    x = 0.5 * lb + 0.5 * ub;

    SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x not in bounds" );

    DenseVector<double> residual( mA * x - b );

    Scalar resNorm = residual.l2Norm();
  
    double tau = 1.0 / residual.l2Norm().getValue<double>();

    SCAI_LOG_INFO ( logger, "start solve, initial resNorm = " << resNorm )

    double gap = 0.0;
    double dualObj = 1.0;

    // Line search parameters

    const double alpha = 0.01;
    const double beta  = 0.5;

    // Parameters for the t-update rule

    const double mu = 2;
    const double s_min = 0.5;

    dualityGap( gap, dualObj, b, x, lb, ub );

    // diagATA = sum(A .* A)';

    mA.reduce( mDiagATA, 1, common::BinaryOp::ADD, common::unary::SQR );

    SCAI_ASSERT_EQ_ERROR( mDiagATA.size(), n, "serious mismatch" )

    SCAI_LOG_INFO( logger, "built diagATA, l2norm = " << mDiagATA.l2Norm()  )

    DenseVector<double> dx( mA.getColDistributionPtr() );

    for ( IndexType iter = 0; iter < mMaxIter; ++iter )
    {
        dx = 0;  // start solution, actually there is no good one

        SCAI_LOG_INFO( logger, "Iter " << iter << ", tau = " << tau << ", gap = " << gap )

        computeSearchDirection( dx, b, x, tau, lb, ub, gap );

        double s = stepSize( b, x, tau, lb, ub, dx, alpha, beta );

        x = x + s * dx;

        SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x no more in boundaries" )

        dualityGap( gap, dualObj, b, x, lb, ub );

        SCAI_LOG_INFO( logger, "gap = " << gap << ", dualObj = " << dualObj )

        if ( gap / common::Math::abs( dualObj ) <= mTolerance )
        {
            // Convergence achieved!
            return;
        }

        if ( s >= s_min )
        {
           tau = common::Math::max( mu * common::Math::min( 2 * n / gap, tau), tau );
        }
    }
}

} // namespace solver

} // namespace scai

