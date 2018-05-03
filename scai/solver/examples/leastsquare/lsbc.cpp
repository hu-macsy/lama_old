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

using namespace scai;
using namespace lama;
using namespace solver;

// return: lb < x < ub, elementwise

typedef DefaultReal ValueType;

bool isInterior( const DenseVector<ValueType>& x,
                 const DenseVector<ValueType>& lb,
                 const DenseVector<ValueType>& ub )
{
    return x.all( common::CompareOp::LT, ub ) && x.all( common::CompareOp::GT, lb );
}

void dualityGap( ValueType& gap, ValueType& dualObj, 
                 const CSRSparseMatrix<ValueType>& A, 
                 const DenseVector<ValueType>& b, 
                 const DenseVector<ValueType>& x, 
                 const DenseVector<ValueType>& lb, 
                 const DenseVector<ValueType>& ub )

{
    // Given x*, compute an estimate for the duality gap by implictly and
    // analytically computing a dual feasible point which converges to the 
    // optimum as x* converges to the optimum for the primal problem

    DenseVector<ValueType> x_s( x - lb );
    DenseVector<ValueType> b_s( b - A * lb );
    DenseVector<ValueType> u_s( ub - lb  );
    DenseVector<ValueType> res( A * x_s - b_s );
    DenseVector<ValueType> kappa( 2 * res );
    DenseVector<ValueType> tmp( kappa * A * Scalar( -1 ) );
    DenseVector<ValueType> mu = tmp;
    mu.setScalar( Scalar( 0 ), common::BinaryOp::MAX );

    Scalar sDualObj = - kappa.dotProduct( kappa ) * 0.25  - kappa.dotProduct( b_s ) - mu.dotProduct( u_s );
    Scalar sGap = res.dotProduct( res ) - sDualObj;

    dualObj = sDualObj.getValue<ValueType>();
    gap = sGap.getValue<ValueType>();
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

        DenseVector<ValueType> tmp( mA * x );
        result += ( 2 * mTau * alpha ) * tmp * mA;
    }

    
    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::ScalarType getValueType() const
    {
        return mA.getValueType();
    }   

private:

    const Matrix& mA;
    const Vector& mD;
    const Scalar mTau;
};

ValueType centralPathObjective(
    const CSRSparseMatrix<ValueType>& A,
    const DenseVector<ValueType>& b,
    const DenseVector<ValueType>& x,
    const ValueType tau,
    const DenseVector<ValueType>& lb,
    const DenseVector<ValueType>& ub )
{
    DenseVector<ValueType> tmp ( x - lb );
    tmp.log(); 
    ValueType s1 = tmp.sum().getValue<ValueType>();
    tmp = ub - x;
    tmp.log(); 
    ValueType s2 = tmp.sum().getValue<ValueType>();
    ValueType barrier = -s1 - s2;
    DenseVector<ValueType> res( A * x - b );
    ValueType dp = res.dotProduct( res ).getValue<ValueType>();
    ValueType value = tau * dp + barrier;

    std::cout << "central path, t = " << tau << ", resnorm = " << common::Math::sqrt( dp ) 
              << ", barrier = " << barrier << ", value = " << value << std::endl;

    return value;
}

ValueType stepSize( const CSRSparseMatrix<ValueType>& A,
                 const DenseVector<ValueType>& b,
                 const DenseVector<ValueType>& x,
                 const ValueType tau,
                 const DenseVector<ValueType>& lb,
                 const DenseVector<ValueType>& ub,
                 const DenseVector<ValueType>& dx,
                 const ValueType alpha,
                 const ValueType beta )
{

//      function [ s ] = step_size (A, b, Vector x, t, l, u, Vector dx, alpha, beta)

// % Gradient g = g(x, t)
// d = 1 ./ (u - x) - 1 ./ (x - l);
// g = 2 * t * A' * (A * x - b) + d;

    DenseVector<ValueType> d1( ub - x ); d1.invert();
    DenseVector<ValueType> d2( x - lb ); d2.invert();
    DenseVector<ValueType> d( d1 - d2 );

    DenseVector<ValueType> g( A * x - b );
    g = 2 * tau * g * A + d;

    IndexType k = 0;

    ValueType objValueAtX = centralPathObjective( A, b, x, tau, lb, ub );

    Scalar objDiff = alpha * g.dotProduct( dx );

    std::cout << "objValue@X = " << objValueAtX << ", diff = " << objDiff << std::endl;

    ValueType s = 1.0;

    while ( true )
    {
        DenseVector<ValueType> x_new( x + s * dx );

        if ( isInterior( x_new, lb, ub ) )
        {
            ValueType objValueAtXNew = centralPathObjective(A, b, x_new, tau, lb, ub );

            ValueType objDiff = objValueAtXNew - objValueAtX;

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
     DenseVector<ValueType>& dx,
     const CSRSparseMatrix<ValueType>& A,
     const DenseVector<ValueType>& b,
     const DenseVector<ValueType>& x,
     const ValueType tau,
     const DenseVector<ValueType>& lb,
     const DenseVector<ValueType>& ub,
     const ValueType gap,
     const DenseVector<ValueType>& diagATA )
{
    // Gradient g = g(x, t)

    // d = 1 ./ (u - x) - 1 ./ (x - l);

    DenseVector<ValueType> d1( ub - x );  
    d1.invert();
    DenseVector<ValueType> d2( x - lb ); 
    d2.invert();
    DenseVector<ValueType> d( d1 - d2 );

    // g = 2 * t * A' * (A * x - b) + d;

    DenseVector<ValueType> tmp( A * x - b );
    DenseVector<ValueType> g( 2 * tau * tmp * A + d );
    g *= -1;

    // Hessian H = H(x, t)
    // D = 1 ./ ((x - l).^2) + 1 ./ ((u - x).^2);
    // H = @(x) 2 * t * A' * (A * x) + D .* x;

    d1 *= d1;
    d2 *= d2;
    DenseVector<ValueType> D( d1 + d2 );

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

    ValueType eps = 1e-14;
    ValueType e = 0.01;

    ValueType normG = g.l2Norm().getValue<ValueType>();
   
    ValueType tol = common::Math::max( eps, common::Math::min( ValueType( 0.1 ), e * gap / normG ) );

    std::cout << "Use relative tol = " << tol << " for CG" << std::endl;

    CriterionPtr criterion1( new IterationCount( 1000 ) );
    NormPtr norm( Norm::create( "L2" ) );   // Norm from factory
    CriterionPtr criterion2( new ResidualThreshold( norm, tol, ResidualThreshold::Relative ) );
    CriterionPtr criterion( new Criterion( criterion1, criterion2, Criterion::OR ) );

    // define common logger to print convergence history

    LoggerPtr logger( new CommonLogger ( "", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly) );

    DenseVector<ValueType> diagonal( 2 * tau * diagATA + D );
    CSRSparseMatrix<ValueType> diagonalMatrix;
    diagonalMatrix.setIdentity( diagonal.getDistributionPtr() );
    diagonalMatrix.setDiagonal( diagonal );

    std::shared_ptr<Jacobi> preconditioner( new Jacobi( "JacobiPreconditioner" ) );
    preconditioner->initialize( diagonalMatrix );

    // Do it with CG

    CG solver( "searchDirectionSolver" );
    solver.setLogger( logger );
    solver.setStoppingCriterion( criterion );
    solver.setPreconditioner( preconditioner );
    solver.initialize( H );
    solver.solve( dx, g );
}

void lsqBox(
    DenseVector<ValueType>& x,
    const CSRSparseMatrix<ValueType>& A,
    const DenseVector<ValueType>& b,
    const DenseVector<ValueType>& lb,
    const DenseVector<ValueType>& ub,
    const ValueType tolerance )
{
    SCAI_REGION ( "lsqBC" )

    const IndexType n = A.getNumColumns();

    SCAI_ASSERT_EQ_ERROR( n, lb.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, ub.size(), "size mismatch" )

    SCAI_ASSERT_EQ_ERROR( b.size(), A.getNumRows(), "size mismatch" )

    SCAI_ASSERT_GT_ERROR( tolerance, ValueType( 0 ), " must be positive" )

    x = 0.5 * lb + 0.5 * ub;

    SCAI_ASSERT_ERROR( isInterior( x, lb, ub ), "x not in bounds" );

    std::cout << "All okay, start least square with box constraints" << std::endl;

    DenseVector<ValueType> residual( A * x - b );

    ValueType tau = 1.0 / residual.l2Norm().getValue<ValueType>();

    ValueType gap = 0.0;
    ValueType dualObj = 1.0;

    // Line search parameters

    const ValueType alpha = 0.01;
    const ValueType beta  = 0.5;

    // Parameters for the t-update rule

    const ValueType mu = 2;
    const ValueType s_min = 0.5;

    const IndexType maxIter = 10000;

    dualityGap( gap, dualObj, A, b, x, lb, ub );

    std::cout << "gap = " << gap << std::endl;

    // diagATA = sum(A .* A)';

    DenseVector<ValueType> diagATA;

    std::cout << "build diagATA" << std::endl;

    A.reduce( diagATA, 1, common::BinaryOp::ADD, common::UnaryOp::SQR );

    std::cout << "diagATA = " << diagATA << std::endl;

    SCAI_ASSERT_EQ_ERROR( diagATA.size(), n, "serious mismatch" )

    // SparseVector<ValueType> col;

    // for( IndexType i = 0; i < n; ++i )
    // {
        // A.getColumn ( col, i );
        // diagATA[i] = col.dotProduct( col ).getValue<ValueType>();
    // }

    std::cout << "norm ATA = " << diagATA.l2Norm() << std::endl;

    DenseVector<ValueType> dx( A.getColDistributionPtr() );

    for ( IndexType iter = 0; iter < maxIter; ++iter )
    {
        dx = 0;  // start solution, actually there is no good one

        std::cout << "Iter " << iter << ", tau = " << tau << ", gap = " << gap << std::endl;

        computeSearchDirection( dx, A, b, x, tau, lb, ub, gap, diagATA);

        ValueType s = stepSize( A, b, x, tau, lb, ub, dx, alpha, beta );

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

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.driver" )
    
    common::Settings::parseArgs( argc, argv );

    CSRSparseMatrix<ValueType> A( "A.mat" );
    DenseVector<ValueType> b ( "b.mat" );
    DenseVector<ValueType> lb ( "lb.mat" );
    DenseVector<ValueType> ub ( "ub.mat" );

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

    // take context as specified by SCAI_CONTEXT

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    A.setContextPtr( ctx );
    b.setContextPtr( ctx );
    ub.setContextPtr( ctx );
    lb.setContextPtr( ctx );

    DenseVector<ValueType> x( ctx );

    ValueType tolerance = 0.01;

    try 
    {
        lsqBox( x, A, b, lb, ub, tolerance );
    }
    catch ( common::Exception& ex )
    {
        std::cout << "Caught exception: " << ex.what() << std::endl;
        std::cout << "Stop execution." << std::endl;
        return 1;
    }

    DenseVector<ValueType> residual( A * x - b );

    std::cout << "res norm = " << residual.l2Norm() << std::endl;

    x.writeToFile( "x.mtx" );
}
