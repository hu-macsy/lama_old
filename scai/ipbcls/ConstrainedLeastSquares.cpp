/**
 * @file ConstrainedLeastSquares.cpp
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
 * @brief Implementation of methods for class to solve LeastSquares with constraints
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#include <scai/hmemo/Context.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include <scai/common/macros/instantiate.hpp>

#include <scai/tracing.hpp>

#include <iostream>

#include <scai/ipbcls/ConstrainedLeastSquares.hpp>
#include <scai/ipbcls/NewtonStepCG.hpp>
#include <scai/ipbcls/CentralPathHessian.hpp>

using std::unique_ptr;
using std::shared_ptr;

namespace scai
{

using common::BinaryOp;
using common::CompareOp;
using common::UnaryOp;
using dmemo::DistributionPtr;
using dmemo::Distribution;
using lama::Vector;
using lama::VectorKind;
using lama::Matrix;
using lama::CentralPathHessian;
using lama::CSRSparseMatrix;
using lama::NormPtr;
using lama::Norm;
using common::Math;
using solver::CriterionPtr;
using solver::CG;
using solver::CommonLogger;
using solver::LogLevel;
using solver::LoggerWriteBehaviour;
using solver::Jacobi;
using solver::Criterion;
using solver::CriterionPtr;
using solver::ResidualThreshold;
using solver::ResidualCheck;
using solver::BooleanOp;
using solver::IterationCount;

namespace ipbcls
{

template<typename ValueType>
bool isInterior( const Vector<ValueType>& x,
                 const Vector<ValueType>& lb,
                 const Vector<ValueType>& ub )
{
    // workaround for: all( x < ub ) && all( x > lb )
    return x.all( CompareOp::LT, ub ) && x.all( CompareOp::GT, lb );
}

template<typename ValueType>
struct DualObjectiveComputation
{
    const Matrix<ValueType>& A;
    const Vector<ValueType>& b;
    const Vector<ValueType>& l;
    const Vector<ValueType>& u;

    unique_ptr<Vector<ValueType>> lambdaPtr;
    unique_ptr<Vector<ValueType>> muPtr;

    DualObjectiveComputation(
        const Matrix<ValueType>& A,
        const Vector<ValueType>& b,
        const Vector<ValueType>& l,
        const Vector<ValueType>& u )

        :   A( A ), b( b ), l( l ), u( u )
    {
        lambdaPtr.reset( A.newSourceVector() );
        muPtr.reset( A.newSourceVector() );
    }

    ValueType operator() ( const Vector<ValueType> & r )
    {
        SCAI_REGION( "ConstrainedLeastSquares.dualityGap" )

        Vector<ValueType>& lambda = *lambdaPtr;
        Vector<ValueType>& mu = *muPtr;

        // Given x*, compute an estimate for the duality gap by implictly and
        // analytically computing a dual feasible point which converges to the
        // optimum as x* converges to the optimum for the primal problem
        // Note that we only need to know r* = A * x* - b.

        // lambda = max(0, 2 A^T r),
        // mu = max(0, - 2 A^T r)
        lambda = 2.0 * transpose( A ) * r;
        mu = -1 * lambda;
        lambda = max( lambda, 0 );
        mu = max ( mu, 0 );

        ValueType rTr = r.dotProduct( r );
        return - rTr - 2.0 * r.dotProduct( b ) + lambda.dotProduct( l ) - mu.dotProduct( u );
    }
};

template<typename ValueType>
ValueType centralPathObjective (
    Vector<ValueType>& tmp,
    const Vector<ValueType>& x,
    const Vector<ValueType>& lb,
    const Vector<ValueType>& ub,
    const ValueType tau,
    const ValueType rTr )
{
    SCAI_REGION( "ConstrainedLeastSquares.centralPathObjective" )
    tmp = x - lb;
    tmp = log( tmp );  // tmp.unaryOpInPlace( common::UNARYOP::LOG );
    const ValueType s1 = tmp.sum();
    tmp = ub - x;
    tmp = log( tmp );  
    const ValueType s2 = tmp.sum();
    const ValueType barrier = - s1 - s2;
    return  tau * rTr + barrier;
}

template<typename ValueType>
struct SearchDirectionComputation
{
    const Matrix<ValueType>& A;
    const Vector<ValueType>& l;
    const Vector<ValueType>& u;

    // diagonal of transpose( A ) * A computed only once by constructor
    std::unique_ptr<const Vector<ValueType>> diagATAPtr;

    // temporary vectors, same size/dist as solution vector
    std::unique_ptr<Vector<ValueType>> gMinusPtr;
    std::unique_ptr<Vector<ValueType>> pPtr;

    CentralPathHessian<ValueType> hessian;
    CSRSparseMatrix<ValueType> pInv;  // diagonal matrix used for preconditioning

    InnerSolverType solverType;

    SearchDirectionComputation( const Matrix<ValueType>& A,
                                const Vector<ValueType>& l,
                                const Vector<ValueType>& u,
                                InnerSolverType solverType )

        : A( A ), l( l ), u( u ), hessian( A ), solverType( solverType )
    {
        const scai::dmemo::DistributionPtr colDist = A.getColDistributionPtr();

        gMinusPtr.reset( A.newSourceVector() );
        pPtr.reset( A.newSourceVector() );

        pInv.setContextPtr( A.getContextPtr() );
        pInv.setIdentity( colDist );

        // diagATA = sum(A .* A)';
        unique_ptr<Vector<ValueType>> diagATAMutablePtr ( A.newSourceVector() );
        A.reduce( *diagATAMutablePtr, 1, BinaryOp::ADD, UnaryOp::SQR );
        diagATAPtr.reset( diagATAMutablePtr.release() );
    }

    // Returns number of iterations
    IndexType operator() (
        Vector<ValueType>& dx,
        const Vector<ValueType> & g,
        const Vector<ValueType> & D,
        const ValueType tau,
        const ValueType gap )
    {
        SCAI_REGION ( "ConstrainedLeastSquares.computeSearchDirection" )

        const Vector<ValueType>& diagATA = *diagATAPtr;
        Vector<ValueType>& gMinus  = *gMinusPtr;
        Vector<ValueType>& p = *pPtr;

        // Need -g as right-hand side for solver
        gMinus = -1 * g;

        // H = 2 * tau * A^T A + D
        hessian.update( D, tau );

        // Set up (diagonal) preconditioner P_inv = 1 ./ (2 * t * diagATA + D)
        // Reuse storage of d to set up preconditioner
        p = 2 * tau * diagATA + D;
        pInv.setDiagonal( p );

        // solver needs a shared pointer for the preconditioner
        shared_ptr<Jacobi<ValueType> > preconditioner( new Jacobi<ValueType>( "JacobiPreconditioner" ) );
        preconditioner->initialize( pInv );

        // TODO: Come up with a better approach to select solver. Currently
        // we allocate a new Solver (and its implicit, internal workspace)
        // for every step computation. We'd naturally like to reuse this,
        // but we need a good way to switch between solvers.
        if ( solverType == InnerSolverType::NewtonStepCG )
        {
            NewtonStepCG<ValueType> solver( &hessian );
            solver.setPreconditioner( preconditioner );

            return solver.computeStep( dx, g );
        }
        else if ( solverType == InnerSolverType::StandardCG )
        {
            // Solve system H dx = -g
            const ValueType eps = 1e-14;
            const ValueType e = 0.01;
            const ValueType normG = g.l2Norm();

            // Note that together with the absolute criterion for the solver,
            // the following constitutes a relative norm norm(r) <= tolerance * norm(g)
            const ValueType tol = Math::max( eps, Math::min( ValueType( 0.1 ), e * gap / normG ) ) * normG;

            CriterionPtr<ValueType> criterion1( new IterationCount<ValueType>( 2000 ) );
            NormPtr<ValueType> norm( Norm<ValueType>::create( "L2" ) );   // Norm from factory
            CriterionPtr<ValueType> criterion2( new ResidualThreshold<ValueType>( norm, tol, ResidualCheck::Absolute ) );
            CriterionPtr<ValueType> criterion( new Criterion<ValueType>( criterion1, criterion2, BooleanOp::OR ) );

            // Allocate a common logger that prints convergenceHistory
            auto logger = std::make_shared<CommonLogger>( "CGLogger: ", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly );

            CG<ValueType> solver( "searchDirectionSolver", logger );
            solver.setStoppingCriterion( criterion );
            solver.setPreconditioner( preconditioner );
            solver.initialize( hessian );
            solver.solve( dx, gMinus );
            return solver.getIterationCount();
        }
        else
        {
            throw std::runtime_error( "Unknown inner solver type requested." );
        }
    }
};

template<typename ValueType>
struct CentralPathComputation
{
    const Matrix<ValueType>& A;
    const Vector<ValueType>& l;
    const Vector<ValueType>& u;

    // tempoarary vectors ( will have same type and same distribution as l, u )

    unique_ptr<Vector<ValueType>> duPtr;
    unique_ptr<Vector<ValueType>> dlPtr;
    unique_ptr<Vector<ValueType>> dPtr;

    CentralPathComputation( const Matrix<ValueType>& A, const Vector<ValueType>& l, const Vector<ValueType>& u )
        :   A( A ), l( l ), u( u )
    {
        duPtr.reset( l.newVector() );
        dlPtr.reset( l.newVector() );
        dPtr.reset( l.newVector() );
    }

    void operator () ( Vector<ValueType>& g, Vector<ValueType>& D, const Vector<ValueType>& r, const Vector<ValueType>& x, const ValueType tau )
    {
        Vector<ValueType>& du = *duPtr;
        Vector<ValueType>& dl = *dlPtr;
        Vector<ValueType>& d  = *dPtr;

        // d = du - dl,
        // where du = 1 ./ (u - x),
        // and   dl = 1 ./ (x - l);
        du = u - x;
        du = 1 / du; 
        dl = x - l;
        dl = 1 / dl;  // not supported yet: dl = 1 / ( x - l );
        d = du - dl;

        // D = du^2 + d^2
        du *= du;
        dl *= dl;
        D = du + dl;

        // Gradient g = g(x, t) : g = 2 * t * A' * (A * x - b) + d;
        g =  2 * tau * transpose( A ) * r + d;
    }
};

template<typename ValueType>
struct StepSizeComputation
{
    const Matrix<ValueType>& A;
    const Vector<ValueType>& l;
    const Vector<ValueType>& u;
    const ValueType alpha;
    const ValueType beta;

    unique_ptr<Vector<ValueType>> xNewPtr;
    unique_ptr<Vector<ValueType>> AdxPtr;
    unique_ptr<Vector<ValueType>> tmpPtr;
    unique_ptr<Vector<ValueType>> rsPtr;

    StepSizeComputation(
        const Matrix<ValueType> & A,
        const Vector<ValueType> & b,
        const Vector<ValueType> & l,
        const Vector<ValueType> & u,
        const ValueType alpha,
        const ValueType beta )

        :   A( A ), l( l ), u( u ), alpha( alpha ), beta( beta )
    {
        xNewPtr.reset ( l.newVector() );
        tmpPtr.reset ( l.newVector() );
        AdxPtr.reset( b.newVector() );
        rsPtr.reset( b.newVector() );
    }

    ValueType operator () (
        const Vector<ValueType>& x,
        const Vector<ValueType>& dx,
        const Vector<ValueType>& r,
        const Vector<ValueType>& g,
        const ValueType tau )
    {
        SCAI_REGION( "ConstrainedLeastSquares.stepSize" )

        Vector<ValueType>& xNew = *xNewPtr;
        Vector<ValueType>& tmp = *tmpPtr;
        Vector<ValueType>& Adx = *AdxPtr;

        // The following identity lets us avoid computing redundant
        // matrix-vector products when evaluating the central path objective
        // for different s:
        // r(s) := A * (x + s * dx) - b = A * x - b + s * dx = r(0) + s * dx = r + s * dx

        Vector<ValueType>& rs = *rsPtr;

        Adx = A * dx;

        ValueType rTr = r.dotProduct( r );
        const ValueType objValueAtX = centralPathObjective( tmp, x, l, u, tau, rTr );
        const ValueType allowedObjDifference = alpha * g.dotProduct( dx );
        ValueType s = 1;

        while ( true )
        {
            xNew = x + s * dx;
            rs = r + s * Adx;
            rTr = rs.dotProduct( rs );

            if ( isInterior( xNew, l, u ) )
            {
                ValueType objValueAtXNew = centralPathObjective( tmp, xNew, l, u, tau, rTr );

                if ( objValueAtXNew <= objValueAtX + s * allowedObjDifference )
                {
                    return s;
                }
            }

            s *= beta;   // make s smaller
        }
    }
};

template<typename ValueType>
bool hasConverged( const ValueType objTolerance,
                   const ValueType resTolerance,
                   const ValueType matrixNormLowerBound,
                   const ValueType bNorm,
                   const ValueType rNorm,
                   const ValueType xNorm,
                   const ValueType gap,
                   const ValueType dualObj )
{
    return dualObj > 0
           ? gap <= objTolerance * dualObj
           : rNorm <= resTolerance * ( matrixNormLowerBound * xNorm + bNorm );
}

template<typename ValueType>
ValueType initialCentralPathParameter( 
    const Vector<ValueType> & x, 
    const Vector<ValueType> & lb, 
    const Vector<ValueType> & ub, 
    const ValueType rTr )
{
    const unique_ptr<Vector<ValueType>> tmp( x.newVector() );

    // The choice of the initial central path parameter is important for the performance of the algorithm:
    // 1. If chosen way too small, one might run into numerical problems due to the requirement of sufficient decrease
    // 2. If chosen too big, the algorithm will convergence extremely slowly because the starting point is very far away
    //    from the central path solution at time t.
    // Recall that the objective function of the central path problem is given by
    //  f(x, t) = t * rTr + phi(x),
    // where r = r(x) and phi(x) = - sum_i [ log(u_i - x_i) ] - sum_i [ log(x_i - l_i) ].
    // We want to pick a t0 such that the minimizer x*(t0) is very close to x*(t).
    // A decent heuristic seems to be
    //   t0 := C * |phi(x*(0))| / dot(r, r),
    // where x*(0) is the minimizer of f(x, 0) and C << 1. Then
    //   f(x*(t0), t0) = C * |phi(x*(0))| + phi(x*(t0)).
    // Since the first term is hopefully small compared to the second,
    // we can hope that x*(t0) is close to x*(0).

    const ValueType initialCentralPathObjAbs = Math::abs( centralPathObjective( *tmp, x, lb, ub, ValueType( 0 ), rTr ) );

    const ValueType c = 0.02;

    // If rTr == 0, it means that the initial residual is zero, which means that
    // the solution has already been found. Hence, we need only make  sure
    // we don't generate an error by dividing by zero here.
    return rTr > 0
           ? c * initialCentralPathObjAbs / rTr
           : c;
}

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ConstrainedLeastSquares<ValueType>::logger, "ConstrainedLeastSquares" )

/* --------------------------------------------------------------------------- */
/*   Constructor                                                               */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
ConstrainedLeastSquares<ValueType>::ConstrainedLeastSquares( const Matrix<ValueType>& A ) :

    mA( A ),
    mObjTolerance( 0.01 ),
    mResTolerance( 1e-6 ),
    mMatrixNormLowerBound ( 0 ),
    mMaxIter( 1000 ),
    mInnerSolverType ( InnerSolverType::NewtonStepCG )
{
    SCAI_LOG_INFO( logger, "Constructed solver for constrained least squares, A = " << A )
}

/* --------------------------------------------------------------------------- */
/*   Setters                                                                   */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ConstrainedLeastSquares<ValueType>::setObjectiveTolerance( ValueType tolerance )
{
    SCAI_ASSERT_GT_ERROR( tolerance, ValueType( 0 ), "tolerance must be positive" )
    mObjTolerance = tolerance;
}

template<typename ValueType>
void ConstrainedLeastSquares<ValueType>::setResidualTolerance( const ValueType tolerance, const ValueType matrixNormLowerBound )
{
    SCAI_ASSERT_GT_ERROR( tolerance, ValueType( 0 ), "tolerance must be positive" )
    SCAI_ASSERT_GE_ERROR( matrixNormLowerBound, ValueType( 0 ), "matrix norm lower bound must be non-negative" )
    mResTolerance = tolerance;
    mMatrixNormLowerBound = matrixNormLowerBound;
}

template<typename ValueType>
ValueType ConstrainedLeastSquares<ValueType>::getObjectiveTolerance() const
{
    return mObjTolerance;
}

template<typename ValueType>
ValueType ConstrainedLeastSquares<ValueType>::getResidualTolerance() const
{
    return mResTolerance;
}

template<typename ValueType>
void ConstrainedLeastSquares<ValueType>::setMaxIter( IndexType maxIter )
{
    mMaxIter = maxIter;
}

template<typename ValueType>
void ConstrainedLeastSquares<ValueType>::setInnerSolverType ( InnerSolverType type )
{
    mInnerSolverType = type;
}

/* --------------------------------------------------------------------------- */
/*   solve                                                                     */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ConstrainedLeastSquares<ValueType>::solve(
    Vector<ValueType>& x,
    const Vector<ValueType>& b,
    const Vector<ValueType>& lb,
    const Vector<ValueType>& ub ) const
{
    SCAI_REGION ( "ConstrainedLeastSquares.solve" )

    const Matrix<ValueType>& A = mA;

    const Distribution& rowDist = A.getRowDistribution();
    const Distribution& colDist = A.getColDistribution();

    SCAI_ASSERT_EQ_ERROR( rowDist, b.getDistribution(), 
                          "Right-hand side b and matrix A are not dimensionally compatible" );
    SCAI_ASSERT_EQ_ERROR( colDist, x.getDistribution(), 
                          "Solution vector x and matrix A are not dimensionally compatible" );
    SCAI_ASSERT_EQ_ERROR( colDist, lb.getDistribution(), 
                          "Lower bound is not dimensionally compatible with solution vector x" )
    SCAI_ASSERT_EQ_ERROR( colDist, ub.getDistribution(), 
                          "Upper bound is not dimensionally compatible" )

    x = 0.5 * lb + 0.5 * ub;

    SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x not in bounds" );

    // Vector<ValueType> is abstract, so we use virtual function newVector to create temporaries

    std::unique_ptr<Vector<ValueType>> dxPtr ( A.newSourceVector() );
    std::unique_ptr<Vector<ValueType>> gPtr ( A.newSourceVector() );
    std::unique_ptr<Vector<ValueType>> DPtr ( A.newSourceVector() );
    std::unique_ptr<Vector<ValueType>> rPtr ( A.newTargetVector() );

    Vector<ValueType>& dx = *dxPtr;
    Vector<ValueType>& g = *gPtr;
    Vector<ValueType>& D = *DPtr;
    Vector<ValueType>& r = *rPtr;

    r =  mA * x - b;
    dx = 0;

    const ValueType rTr_initial = r.dotProduct( r );
    const ValueType bNorm = b.l2Norm();
    ValueType tau = initialCentralPathParameter( x, lb, ub, rTr_initial );
    IndexType totalInnerIterations = 0;

    // Line search parameters
    const ValueType alpha = 0.01;
    const ValueType beta  = 0.5;

    // Parameters for the tau-update rule
    const ValueType mu = 2;
    const ValueType s_min = 0.5;

    DualObjectiveComputation<ValueType> dualObjective( A, b, lb, ub );
    CentralPathComputation<ValueType> computeCentralPathComponents( A, lb, ub );
    StepSizeComputation<ValueType> stepSize( A, b, lb, ub, alpha, beta );
    SearchDirectionComputation<ValueType> computeSearchDirection( A, lb, ub, mInnerSolverType );

    ValueType bestDualObj = -1e300;    // is this okay for ValueType == float

    for ( IndexType iter = 0; iter < mMaxIter; ++iter )
    {
        const ValueType rTr = r.dotProduct( r );
        const ValueType rNorm = Math::sqrt( rTr );
        const ValueType dualObj = dualObjective( r );
        bestDualObj = Math::max( dualObj, bestDualObj );
        const ValueType gap = rTr - bestDualObj;

        SCAI_LOG_INFO( logger, "Iteration " << iter << ": tau = " << tau << ", rTr = " << rTr << ", gap = " << gap )

        if ( hasConverged( getObjectiveTolerance(), getResidualTolerance(), mMatrixNormLowerBound, bNorm,
                           rNorm, x.l2Norm(), gap, dualObj ) )
        {
            SCAI_LOG_INFO( logger, "Tolerance achieved, stop at iteration = " << iter
                           << ". Total inner iterations: " << totalInnerIterations )
            return;
        }

        computeCentralPathComponents( g, D, r, x, tau );
        SCAI_LOG_INFO( logger, "computeSearchDirection" )
        const IndexType innerIterations = computeSearchDirection( dx, g, D, tau, gap );
        SCAI_LOG_INFO( logger, "computeSearchDirection, #iter = " << innerIterations )
        const ValueType s = stepSize( x, dx, r, g, tau );

        x = x + s * dx;
        r = A * x - b;

        if ( s >= s_min )
        {
            dx = 0;
            const IndexType n = colDist.getGlobalSize();
            tau = Math::max( mu * Math::min( ValueType( 2 * n ) / gap, tau ), tau );
        }
        else
        {
            dx = s * dx;
        }

        totalInnerIterations += innerIterations;
        SCAI_LOG_DEBUG( logger, "Inner iteration count = " << innerIterations )
        SCAI_LOG_DEBUG( logger, "Step size s = " << s )
        SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "Violation of invariant: x is not in the interior" )
    }

    SCAI_LOG_INFO( logger, "Solver did not converge within maxmimum number of iterations."
                   << "Max iterations: " << mMaxIter
                   << ". Total inner iterations: " << totalInnerIterations )
}

// instantiation only for real types ( as < / > constraints are not defined for complex )

SCAI_COMMON_INST_CLASS( ConstrainedLeastSquares, SCAI_REAL_TYPES_HOST )

}

}
