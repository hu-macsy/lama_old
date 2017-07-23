/**
 * @file lsbc.cpp
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
 * @brief Example of least square problem with boundary conditions
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 15.11.2016
 */

#include <scai/solver/examples/eigenvalue/AbstractMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/common/Math.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/MINRES.hpp>
#include <scai/solver/GMRES.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>

#include "DiagonalSolver.hpp"

#include <iostream>

using namespace scai;
using namespace lama;
using namespace solver;

// return: lb < x < ub, elementwise

bool isInterior( const DenseVector<double>& x,
                 const DenseVector<double>& lb,
                 const DenseVector<double>& ub )
{
    // SCAI_ASSERT_EQ_ERROR( x.getDistribution(), lb.getDistribution(), "distributon mismatch" )
    // SCAI_ASSERT_EQ_ERROR( x.getDistribution(), ub.getDistribution(), "distributon mismatch" )

    // bool okay = true;

    // hmemo::ReadAccess<double> rX( x.getLocalValues() );
    // hmemo::ReadAccess<double> rLB( lb.getLocalValues() );
    // hmemo::ReadAccess<double> rUB( ub.getLocalValues() );

    // for ( IndexType i = 0; i < rX.size(); ++i )
    // {
         // okay = okay && ( rLB[i] < rX[i] ) && ( rX[i] < rUB[i] );
    // }

    // okay = x.getDistribution().getCommunicator().all( okay );

    // return okay;

    return x.all( common::binary::LT, ub ) && x.all( common::binary::GT, lb );
}

void dualityGap( double& gap, double& dualObj, 
                 const CSRSparseMatrix<double>& A, 
                 const DenseVector<double>& b, 
                 const DenseVector<double>& x, 
                 const DenseVector<double>& lb, 
                 const DenseVector<double>& ub )

{
    // Given x*, compute an estimate for the duality gap by implictly and
    // analytically computing a dual feasible point which converges to the 
    // optimum as x* converges to the optimum for the primal problem

    DenseVector<double> x_s( x - lb );
    DenseVector<double> b_s( b - A * lb );
    DenseVector<double> u_s( ub - lb  );
    DenseVector<double> res( A * x_s - b_s );
    DenseVector<double> kappa( 2 * res );
    DenseVector<double> tmp( kappa * A * Scalar( -1 ) );
    DenseVector<double> mu = tmp;
    mu.setScalar( Scalar( 0 ), common::binary::MAX );

    Scalar sDualObj = - kappa.dotProduct( kappa ) * 0.25  - kappa.dotProduct( b_s ) - mu.dotProduct( u_s );
    Scalar sGap = res.dotProduct( res ) - sDualObj;

    dualObj = sDualObj.getValue<double>();
    gap = sGap.getValue<double>();
}

class HessianMatrix : public AbstractMatrix 
{

public:

    HessianMatrix ( const Matrix& A, const Vector& D, const Scalar& tau ) :

        AbstractMatrix( A.getColDistributionPtr(), A.getColDistributionPtr() ),

        mA( A ),
        mD( D ),
        mTau( tau )
    {
    }

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_LOG_INFO( logger, "matrixTimesVector, mA = " << mA )

        result = mD * x;
        result += beta * y;

        DenseVector<double> tmp( mA * x );
        result += ( 2 * mTau * alpha ) * tmp * mA;
    }

    
    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::scalar::ScalarType getValueType() const
    {
        return mA.getValueType();
    }   

private:

    const Matrix& mA;
    const Vector& mD;
    const Scalar mTau;
};

double centralPathObjective(
    const CSRSparseMatrix<double>& A,
    const DenseVector<double>& b,
    const DenseVector<double>& x,
    const double tau,
    const DenseVector<double>& lb,
    const DenseVector<double>& ub )
{
    DenseVector<double> tmp ( x - lb );
    tmp.log(); 
    double s1 = tmp.sum().getValue<double>();
    tmp = ub - x;
    tmp.log(); 
    double s2 = tmp.sum().getValue<double>();
    double barrier = -s1 - s2;
    DenseVector<double> res( A * x - b );
    double dp = res.dotProduct( res ).getValue<double>();
    double value = tau * dp + barrier;

    std::cout << "central path, t = " << tau << ", resnorm = " << common::Math::sqrt( dp ) 
              << ", barrier = " << barrier << ", value = " << value << std::endl;

    return value;
}

double stepSize( const CSRSparseMatrix<double>& A,
                 const DenseVector<double>& b,
                 const DenseVector<double>& x,
                 const double tau,
                 const DenseVector<double>& lb,
                 const DenseVector<double>& ub,
                 const DenseVector<double>& dx,
                 const double alpha,
                 const double beta )
{

//      function [ s ] = step_size (A, b, Vector x, t, l, u, Vector dx, alpha, beta)

// % Gradient g = g(x, t)
// d = 1 ./ (u - x) - 1 ./ (x - l);
// g = 2 * t * A' * (A * x - b) + d;

    DenseVector<double> d1( ub - x ); d1.invert();
    DenseVector<double> d2( x - lb ); d2.invert();
    DenseVector<double> d( d1 - d2 );

    DenseVector<double> g( A * x - b );
    g = 2 * tau * g * A + d;

    IndexType k = 0;

    double objValueAtX = centralPathObjective( A, b, x, tau, lb, ub );

    Scalar objDiff = alpha * g.dotProduct( dx );

    std::cout << "objValue@X = " << objValueAtX << ", diff = " << objDiff << std::endl;

    double s = 1.0;

    while ( true )
    {
        DenseVector<double> x_new( x + s * dx );

        if ( isInterior( x_new, lb, ub ) )
        {
            double objValueAtXNew = centralPathObjective(A, b, x_new, tau, lb, ub );

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

void computeSearchDirection(
     DenseVector<double>& dx,
     const CSRSparseMatrix<double>& A,
     const DenseVector<double>& b,
     const DenseVector<double>& x,
     const double tau,
     const DenseVector<double>& lb,
     const DenseVector<double>& ub,
     const double gap,
     const DenseVector<double>& diagATA )
{
    // Gradient g = g(x, t)

    // d = 1 ./ (u - x) - 1 ./ (x - l);

    DenseVector<double> d1( ub - x );  
    d1.invert();
    DenseVector<double> d2( x - lb ); 
    d2.invert();
    DenseVector<double> d( d1 - d2 );

    // g = 2 * t * A' * (A * x - b) + d;

    DenseVector<double> tmp( A * x - b );
    DenseVector<double> g( 2 * tau * tmp * A + d );
    g *= -1;

    // Hessian H = H(x, t)
    // D = 1 ./ ((x - l).^2) + 1 ./ ((u - x).^2);
    // H = @(x) 2 * t * A' * (A * x) + D .* x;

    d1 *= d1;
    d2 *= d2;
    DenseVector<double> D( d1 + d2 );

    HessianMatrix H( A, D, tau );

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

    std::cout << "Use relative tol = " << tol << " for CG" << std::endl;

    CriterionPtr criterion1( new IterationCount( 1000 ) );
    NormPtr norm( Norm::create( "L2" ) );   // Norm from factory
    CriterionPtr criterion2( new ResidualThreshold( norm, tol, ResidualThreshold::Relative ) );
    CriterionPtr criterion( new Criterion( criterion1, criterion2, Criterion::OR ) );

    // define common logger to print convergence history

    LoggerPtr logger( new CommonLogger ( "", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly) );

    DenseVector<double> diagonal( 2 * tau * diagATA + D );
    CSRSparseMatrix<double> diagonalMatrix;
    diagonalMatrix.setIdentity( diagonal.getDistributionPtr() );
    diagonalMatrix.setDiagonal( diagonal );

    common::shared_ptr<Jacobi> preconditioner( new Jacobi( "JacobiPreconditioner" ) );
    preconditioner->initialize( diagonalMatrix );

    // common::shared_ptr<DiagonalSolver> preconditioner( new DiagonalSolver( "JacobiPreconditioner" ) );
    // preconditioner->initialize( diagonal );

    // Do it with CG

    CG solver( "searchDirectionSolver" );
    solver.setLogger( logger );
    solver.setStoppingCriterion( criterion );
    solver.setPreconditioner( preconditioner );
    solver.initialize( H );
    solver.solve( dx, g );
}

void lsqBox(
    DenseVector<double>& x,
    const CSRSparseMatrix<double>& A,
    const DenseVector<double>& b,
    const DenseVector<double>& lb,
    const DenseVector<double>& ub,
    const double tolerance )
{
    const IndexType n = A.getNumColumns();

    SCAI_ASSERT_EQ_ERROR( n, lb.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, ub.size(), "size mismatch" )

    SCAI_ASSERT_EQ_ERROR( b.size(), A.getNumRows(), "size mismatch" )

    SCAI_ASSERT_GT_ERROR( tolerance, double( 0 ), " must be positive" )

    x = 0.5 * lb + 0.5 * ub;

    SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x not in bounds" );

    std::cout << "All okay, start least square with box constraints" << std::endl;

    DenseVector<double> residual( A * x - b );

    double tau = 1.0 / residual.l2Norm().getValue<double>();

    double gap = 0.0;
    double dualObj = 1.0;

    // Line search parameters

    const double alpha = 0.01;
    const double beta  = 0.5;

    // Parameters for the t-update rule

    const double mu = 2;
    const double s_min = 0.5;

    const IndexType maxIter = 10000;

    dualityGap( gap, dualObj, A, b, x, lb, ub );

    std::cout << "gap = " << gap << std::endl;

    // diagATA = sum(A .* A)';

    DenseVector<double> diagATA;

    std::cout << "build diagATA" << std::endl;

    A.reduce( diagATA, 1, common::binary::ADD, common::unary::SQRT );

    std::cout << "diagATA = " << diagATA << std::endl;

    SCAI_ASSERT_EQ_ERROR( diagATA.size(), n, "serious mismatch" )

    // SparseVector<double> col;

    // for( IndexType i = 0; i < n; ++i )
    // {
        // A.getColumn ( col, i );
        // diagATA[i] = col.dotProduct( col ).getValue<double>();
    // }

    std::cout << "norm ATA = " << diagATA.l2Norm() << std::endl;

    DenseVector<double> dx( A.getColDistributionPtr() );

    for ( IndexType iter = 0; iter < maxIter; ++iter )
    {
        dx = 0;  // start solution, actually there is no good one

        std::cout << "Iter " << iter << ", tau = " << tau << ", gap = " << gap << std::endl;

        computeSearchDirection( dx, A, b, x, tau, lb, ub, gap, diagATA);

        double s = stepSize( A, b, x, tau, lb, ub, dx, alpha, beta );

        std::cout << "Step size = " << s << std::endl;

        x = x + s * dx;

        SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x no more in boundaries" )

        dualityGap( gap, dualObj, A, b, x, lb, ub );

        std::cout << "gap = " << gap << ", dualObj = " << dualObj << std::endl;

        // TODO: Is it correct to take abs of dual_obj here?

        if ( gap / common::Math::abs( dualObj ) <= tolerance )
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

int main( int, char** )
{
    CSRSparseMatrix<double> A( "A.mat" );
    DenseVector<double> b ( "b.mat" );
    DenseVector<double> lb ( "lb.mat" );
    DenseVector<double> ub ( "ub.mat" );

    std::cout << "A = " << A << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "lb = " << lb << std::endl;
    std::cout << "ub = " << ub << std::endl;

    if ( true )
    {
        A.writeToFile( "A.mtx" );
        b.writeToFile( "b.mtx" );
        lb.writeToFile( "lb.mtx" );
        ub.writeToFile( "ub.mtx" );
    }

    DenseVector<double> x;

    DenseVector<double> range( ub - lb );

    std::cout << "range = " << range.min() << " - " << range.max() << std::endl;

    double tolerance = 0.001;

    lsqBox( x, A, b, lb, ub, tolerance );

    DenseVector<double> residual( A * x - b );

    std::cout << "res norm = " << residual.l2Norm() << std::endl;
}
