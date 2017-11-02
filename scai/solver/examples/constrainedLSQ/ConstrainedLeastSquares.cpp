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
 * @author Thomas Brandes, Andreas Borgen Longva
 * @date 21.07.2017
 */

#include <scai/hmemo/Context.hpp>

#include <scai/solver/CG.hpp>
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

bool isInterior( const Vector& x,
                 const Vector& lb,
                 const Vector& ub )
{
    return x.all( common::binary::LT, ub ) && x.all( common::binary::GT, lb );
}

void ConstrainedLeastSquares::dualityGap( 

    Scalar& gap, 
    Scalar& dualObj, 
    const Vector& b, 
    const Vector& x, 
    const Vector& lb, 
    const Vector& ub )

{
    SCAI_REGION( "lsqb.dualityGap" )

    // Given x*, compute an estimate for the duality gap by implictly and
    // analytically computing a dual feasible point which converges to the 
    // optimum as x* converges to the optimum for the primal problem

    Vector& d1    = *d1Ptr;
    Vector& d2    = *d2Ptr;

    Vector& kappa = *kappaPtr;

    d1 = x - lb;
    _VectorPtr bsPtr( mA.newVector( 0 ) );
    Vector& b_s = *bsPtr;
    b_s = b - mA * lb;
    d2 = ub - lb;
    Vector& res = *resPtr;
    res = mA * d1 - b_s;
    kappa = 2 * res;

    Vector& m = *mPtr;
    m = ( -1 ) * kappa * mA;

    m.setScalar( Scalar( 0 ), common::binary::MAX );

    dualObj = - kappa.dotProduct( kappa ) * 0.25  - kappa.dotProduct( b_s ) - m.dotProduct( d2 );
    gap     = res.dotProduct( res ) - dualObj;

    SCAI_LOG_INFO( logger, "dualObj = " << dualObj << ", gap = " << gap )
}

Scalar ConstrainedLeastSquares::centralPathObjective (
    const Vector& b,
    const Vector& x,
    const Scalar tau,
    const Vector& lb,
    const Vector& ub )
{
    SCAI_REGION( "lsqb.centralPathObjective" )

    Vector& d1 = *d1Ptr;
    Vector& res = *resPtr;

    d1 = x - lb;  d1.log(); 
    Scalar s1 = d1.sum();

    d1 = ub - x; d1.log(); 
    Scalar s2 = d1.sum();

    Scalar barrier = -s1 - s2;

    res =  mA * x - b;
    Scalar dp = res.dotProduct( res );
    Scalar value = tau * dp + barrier;

    SCAI_LOG_INFO( logger, "central path, t = " << tau << ", resnorm = " << sqrt( dp ) 
                            << ", barrier = " << barrier << ", value = " << value )

    return value;
}

Scalar ConstrainedLeastSquares::stepSize( 
    const Vector& b,
    const Vector& x,
    const Scalar tau,
    const Vector& lb,
    const Vector& ub,
    const Vector& dx,
    const Scalar alpha,
    const Scalar beta )
{
    SCAI_REGION( "lsqb.stepSize" )

//      function [ s ] = step_size (A, b, Vector x, t, l, u, Vector dx, alpha, beta)

// % Gradient g = g(x, t)
// d = 1 ./ (u - x) - 1 ./ (x - l);
// g = 2 * t * A' * (A * x - b) + d;

    Vector& d1 = *d1Ptr;
    Vector& d2 = *d2Ptr;
    Vector& g  = *gPtr;
    Vector& x_new = *xNewPtr;

    Vector& res = *resPtr;   // use temporary here

    d1 = ub - x; d1.invert();
    d2 = x -lb;  d2.invert();

    d1 -= d2;

    res = mA * x - b;

    g = 2 * tau * res * mA + d1;

    IndexType k = 0;

    Scalar objValueAtX = centralPathObjective( b, x, tau, lb, ub );

    Scalar objDiff = alpha * g.dotProduct( dx );

    std::cout << "objValue@X = " << objValueAtX << ", diff = " << objDiff << std::endl;

    Scalar s = 1;

    while ( true )
    {
        x_new = x + s * dx;

        if ( isInterior( x_new, lb, ub ) )
        {
            Scalar objValueAtXNew = centralPathObjective(  b, x_new, tau, lb, ub );

            Scalar objDiff = objValueAtXNew - objValueAtX;

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
     Vector& dx,
     const Vector& b,
     const Vector& x,
     const Scalar tau,
     const Vector& lb,
     const Vector& ub,
     const Scalar gap )
{
    SCAI_REGION ( "clsq.computeSearchDirection" )

    Vector& d1 = *d1Ptr;
    Vector& d2 = *d2Ptr;
    Vector& d  = *dPtr;
    Vector& D  = *DPtr;
    Vector& g  = *gPtr;

    // d = 1 ./ (u - x) - 1 ./ (x - l);

    d1 = ub - x;  
    d1.invert();
    d2 = x - lb; 
    d2.invert();
    d = d1 - d2;

    // Gradient g = g(x, t) : g = 2 * t * A' * (A * x - b) + d;

    Vector& tmp = *resPtr;
    tmp = mA * x - b;

    g =  2 * tau * tmp * mA + d;
    g *= -1;

    // Hessian H = H(x, t)
    // D = 1 ./ ((x - l).^2) + 1 ./ ((u - x).^2);
    // H = @(x) 2 * t * A' * (A * x) + D .* x;

    d1 *= d1;
    d2 *= d2;
    D = d1 + d2;

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

    Scalar eps = 1e-14;
    Scalar e = 0.01;

    Scalar normG = g.l2Norm().getValue<double>();
   
    Scalar tol = max( eps, min( Scalar( 0.1 ), e * gap / normG ) );

    SCAI_LOG_DEBUG( logger, "Use relative tol = " << tol << " for CG" )

    CriterionPtr criterion1( new IterationCount( 2000 ) );
    NormPtr norm( Norm::create( "L2" ) );   // Norm from factory
    CriterionPtr criterion2( new ResidualThreshold( norm, tol, ResidualThreshold::Relative ) );
    CriterionPtr criterion( new Criterion( criterion1, criterion2, Criterion::OR ) );

    // define common logger to print convergence history

    d = 2 * tau * mDiagATA + D;

    mDiagonalMatrix.setDiagonal( d );

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

/* --------------------------------------------------------------------------- */
/*   Constructor                                                               */
/* --------------------------------------------------------------------------- */

ConstrainedLeastSquares::ConstrainedLeastSquares( const lama::_Matrix& A ) :

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

    // allocate runtime vectors

    resPtr.reset( mA.newVector( 0 ) );  // size of residual = #rows of A
    kappaPtr.reset( mA.newVector( 0 ) );

    d1Ptr.reset( mA.newVector( 1 ) );  // size of d1 = size of x = size of #cols of A
    d2Ptr.reset( mA.newVector( 1 ) );  // size of d1 = size of x = size of #cols of A
    dxPtr.reset( mA.newVector( 1 ) );  // deltaX
    dPtr.reset( mA.newVector( 1 ) );   // d
    DPtr.reset( mA.newVector( 1 ) );   // D
    gPtr.reset( mA.newVector( 1 ) );   // g
    mPtr.reset( mA.newVector( 1 ) );  // mu
    xNewPtr.reset( mA.newVector( 1 ) );  // mu
}

/* --------------------------------------------------------------------------- */
/*   Setters                                                                   */
/* --------------------------------------------------------------------------- */

void ConstrainedLeastSquares::setTolerance( Scalar tolerance )
{
    SCAI_ASSERT_GT_ERROR( tolerance, Scalar( 0 ), "tolerance must be positive" )
    mTolerance = tolerance;
}

void ConstrainedLeastSquares::setMaxIter( IndexType maxIter )
{
    mMaxIter = maxIter;
}

/* --------------------------------------------------------------------------- */
/*   solve                                                                     */
/* --------------------------------------------------------------------------- */

void ConstrainedLeastSquares::solve(
    Vector& x,
    const Vector& b,
    const Vector& lb,
    const Vector& ub )
{
    SCAI_REGION ( "clsq.solve" )

    const IndexType n = mA.getNumColumns();

    SCAI_ASSERT_EQ_ERROR( n, lb.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, ub.size(), "size mismatch" )

    SCAI_LOG_INFO( logger, "sizes: b = " << b.size() << "#rows A = " << mA.getNumRows() )
    SCAI_LOG_INFO( logger, "b = " << b )

    SCAI_ASSERT_EQ_ERROR( b.size(), mA.getNumRows(), "size mismatch" )

    x = 0.5 * lb + 0.5 * ub;

    SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x not in bounds" );

    _VectorPtr residualPtr( mA.newVector( 0 ) );   // size = #rows
    *residualPtr = 0;

    SCAI_LOG_INFO( logger, "norm b = " << b.l2Norm() )

    *residualPtr =  mA * x - b;

    Scalar resNorm = residualPtr->l2Norm();

    SCAI_LOG_INFO( logger, "initial res norm = " << resNorm )

    Scalar tau = Scalar( 1 ) / resNorm;

    Scalar gap = 0.0;
    Scalar dualObj = 1.0;

    // Line search parameters

    const Scalar alpha = 0.01;
    const Scalar beta  = 0.5;

    // Parameters for the t-update rule

    const Scalar mu = 2;
    const double s_min = 0.5;

    dualityGap( gap, dualObj, b, x, lb, ub );

    // diagATA = sum(A .* A)';

    mA.reduce( mDiagATA, 1, common::binary::ADD, common::unary::SQR );

    SCAI_ASSERT_EQ_ERROR( mDiagATA.size(), n, "serious mismatch" )

    SCAI_LOG_INFO( logger, "built diagATA, l2norm = " << mDiagATA.l2Norm()  )

    Vector& dx = *dxPtr;

    for ( IndexType iter = 0; iter < mMaxIter; ++iter )
    {
        dx = 0;  // start solution, actually there is no good one

        SCAI_LOG_INFO( logger, "Iter " << iter << ", tau = " << tau << ", gap = " << gap )

        computeSearchDirection( dx, b, x, tau, lb, ub, gap );

        Scalar s = stepSize( b, x, tau, lb, ub, dx, alpha, beta );

        x = x + s * dx;

        SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x no more in boundaries" )

        dualityGap( gap, dualObj, b, x, lb, ub );

        SCAI_LOG_INFO( logger, "gap = " << gap << ", dualObj = " << dualObj )

        if ( gap / abs( dualObj ) < mTolerance )
        {
            // Convergence achieved!
            return;
        }

        if ( s > s_min )
        {
           tau = max( mu * min( Scalar( 2 * n ) / gap, tau), tau );
        }
    }
}

} // namespace solver

} // namespace scai

