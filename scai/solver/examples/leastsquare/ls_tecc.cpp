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
 * @author Thomas Brandes, Andreas Borgen Langva, Dustin Feld, Lauretta Schubert
 * @date 15.11.2016
 */

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/common/Math.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/solver/CG.hpp>

#include <iostream>

using namespace scai::solver;
using namespace scai::lama;

typedef RealType ValueType;

void joinMatrix( CSRSparseMatrix<ValueType>& result, const CSRSparseMatrix<ValueType>& a, const CSRSparseMatrix<ValueType>& b )
{
    SCAI_ASSERT_EQ_ERROR( a.getNumColumns(), b.getNumColumns(), "joined matrices must have same number of columns" );

    typedef std::shared_ptr<scai::lama::_MatrixStorage> StoragePtr;

    StoragePtr shared_ptrA( a.getLocalStorage().copy() );
    StoragePtr shared_ptrB( b.getLocalStorage().copy() );

    std::vector<StoragePtr> bothMatrices;

    bothMatrices.push_back( shared_ptrA );
    bothMatrices.push_back( shared_ptrB );

    scai::lama::CSRStorage<ValueType> joinedStorage;

    joinedStorage.rowCat( bothMatrices );

    result.assign( joinedStorage );
}

void setupSmoothMatrix( CSRSparseMatrix<ValueType>& L, const IndexType yMax, const IndexType zMax )
{
    common::Stencil2D stencil( 5 );

    StencilMatrix<ValueType> stencilMatrix( stencil, common::Grid2D( yMax, zMax ) );

    L = stencilMatrix;
}

void zeroExtend( DenseVector<ValueType>& T_ext, 
                 const DenseVector<ValueType>& T, const IndexType nZeros )
{
    T_ext.allocate(	T.size() + nZeros );

    hmemo::WriteAccess<ValueType> wT( T_ext.getLocal() );
    hmemo::ReadAccess<ValueType> rT( T.getLocal() );

    for ( IndexType i = 0; i < rT.size(); ++i )
    {
        wT[i] = rT[i];
    }
}

void computeSearchDirection( 
     DenseVector<ValueType>& dx,
     const CSRSparseMatrix<ValueType>& A,
     const DenseVector<ValueType>& lb,
     const DenseVector<ValueType>& ub,
     const ValueType tau,
     const DenseVector<ValueType>& ATA_diag,

{
    // Gradient g = g(x, t)

    // d = 1 ./ (u - x) - 1 ./ (x - l);

    DenseVector<ValueType> d1 = ub - x; d1.invert();
    DenseVector<ValueType> d2 = x - lb; d2.invert();
    DenseVector<ValueType> d  = d1 - d2;
    
    // g = 2 * t * A' * (A * x - b) + d;

    DenseVector<ValueType> g( A * x - b );
    g = 2 * tau * g * A + d;

    // % Hessian H = H(x, t)
    // D = 1 ./ ((x - l).^2) + 1 ./ ((u - x).^2);
    // H = @(x) 2 * t * A' * (A * x) + D .* x;

   // % Preconditioner for H (note that it is diagonal)
   // P_inv = 1 ./ (2 * t * ATA_diag + D);

   // % Solve H dx = -g
   // % e is just a multiplier for the tolerance
   // e = 0.01;
   // tol = max(1e-14, min(0.1, e * gap / norm(g)));

   /// relative tolerance, 

   // [dx, ~, ~, ~] = pcg(H, - g, tol, 1000, @(x) P_inv .* x);
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

    DenseVector<ValueType> d1 = ub - x; d1.invert();
    DenseVector<ValueType> d2 = x - lb; d2.invert();
    DenseVector<ValueType> d  = d1 - d2;
    
    DenseVector<ValueType> g( A * x - b );
    g = 2 * tau * g * A + d;

    IndexType k = 0;

    ValueType objValueAtX = centralPathObjective( A, b, x, tau, lb, ub );

    while ( true )
    {
        ValueType s = beta ^ k;

        Vector x_new = x + s * dx;

        ValueType objValueAtXNew = centralPathObjective(A, b, x_new, tau, lb, ub );

        if ( !isInterior( x_new, lb, ub ) )
        {
           
        }
    
        if ( isInterior( ... ) &&  ( objValueAtXNew <= objValueAtX + alpha * s * g.dotProduct( dx) )
        {
            return s;
        }
    
        k = k + 1;
    }
}

/** This method finds the min of || A x - b ||_2 with lb <= x <= ub 
 *
 *  ToDo: how to specifiy lb[i] == -Inf or ub[i] == Inf
 *
 */
void lsqBox( DenseVector<ValueType>& x,
             const CSRSparseMatrix<ValueType>& A, 
             const DenseVector<ValueType>& b, 
             const DenseVector<ValueType>& lb,
             const DenseVector<ValueType>& ub,
             const ValueType tolerance )
{

    SCAI_ASSERT_EQ_ERROR( A.getNumColumns(), lb.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( A.getNumColumns(), ub.size(), "size mismatch" )

    SCAI_ASSERT_EQ_ERROR( b.size(), A.getNumRows(), "size mismatch" )

    // ToDo: lb < ub 

    SCAI_ASSERT_GT_ERROR( tolerance, ValueType( 0 ), " must be positive" )

    const IndexType n = A.getNumColumns();

    // Line search parameters, recommended values

    ValueType alpha = 0.01;
    ValurType beta  = 0.5;

    // Parameters for the t-update rule

    ValueType mu = 2;
    ValueType s_min = 0.5;

    //  Set the initial point to the middle of the interior, which lets us
   //  focus on homing in on the solution rather than moving towards the
   //  interior.

    x = 0.5 * lb + 0.5 * ub;

    // TODO: There may be a more sensible way to initialize t, the interior point

    DenseVector<ValueType> residual( A * x - b );

    Scalar t = 1 / residual.l2Norm();

    dualityGap( gap, dualOjb, A, b, x, lb, ub );

    // ATA_diag = sum(A .* A)';

    DenseVector<ValueType> ATA_diag( n );

    SparseVector<ValueType> col;

    for( IndexType i = 0; i < n; ++i )
    {
        A.getColumn ( col, i );
        ATA_diag[i] = col.dotProdcut( col );
    }

    for ( IndexType iter = 0; iter < MaxIter; ++iter )
    {
        computeSearchDirection( dx, A, b, x, t, l, u, gap, ATA_diag);

        s = step_size( A, b, x, t, l, u, dx, alpha, beta );
        x = x + s * dx;
    
        dualityGap( gap, dualOjb, A, b, x, lb, ub );
    
        assert(all(l < x) & all(u > x));
    
        // TODO: Is it correct to take abs of dual_obj here?

        if ( gap / abs(dual_obj) <= tolerance )
        {
            // Convergence achieved!
            return
        }
    
        if ( s >= s_min )
        {
           t = max(mu * min(2 * n / gap, t), t);
        }
    }
}

class HessianMatrix : AbstractMatrix ( const Matrix& A, const Vector& D, const Scalar& t )
{
    HessianMatrix ( const Matrix& A, const Vector& D, const Scalar& t )
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

        result = elemWiseMult( D , x );
        tmp = A * x;
        result = ( 2 * t ) * tmp * A  + result + D * x;
    }

private

    const Matrix& mA;
    const Vector& mD;
    const Scalar mT;
};

// retrn: bool inBounds

bool inBounds( const DenseVector<ValueType>& x,
               const DenseVector<ValueType>& lb,
               const DenseVector<ValueType>& ub )
{
     bool okay = true;

     ReadAccess<ValueType> rX( x.local() );
     ReadAccess<ValueType> rLB( lb.local() );
     ReadAccess<ValueType> rUB( ub.local() );

     for ( IndexType i = 0; i < rX.size(), ++i )
     {
          okay = okay && ( rLB[i] < rX[i] ) && ( rX[i] < ruB[i] );
     }

     okay = x.getDistribution().getCommunicator().all( okay );

     return okay;
}

ValueTypValueType centralPathObjective( 
    const CSRSparseMatrix<ValueType>& A,
    const DenseVector<ValueType>& b,
    const Scalar tau,
    const DenseVector<ValueType>& lb,
    const DenseVector<ValueType>& ub )
{
    DenseVector<ValueType> tmp = x - lb;
    tmp.log();
    ValueType s1 = tmp.sum();
    tmp = ub - x;
    tmp.log();
    ValueType s2 = tmp.sum();
    ValueType barrier = -s1 - s2;
    DenseVector<ValueType> res( A * x - b );
    value = tau * res.dotProduct( res ) + barrier;
    return value;
}

void dualityGap( Scalar& gap, Scalar& dualObj, const Matrix& A, cont Vector& b, const Vector& x, const Vector& l, const Vector& u )

{
    // Given x*, compute an estimate for the duality gap by implictly and
    // analytically computing a dual feasible point which converges to the 
    // optimum as x* converges to the optimum for the primal problem

    VectorPtr b_s ( b.newVector() );

    *b_s = b - A * l;

    Vector u_s = u - l;

    Vector r( b.newVector() );

    *r  = A * x - b;

    Vector kappa( b.newVector() );

    *kappa = 2 * r;

    Vector tmp = - kappa * A;

    DenseVector<ValueType>  zeros( x.size(), 0 );
    Vector mu = max ( zeros, tmp );
    
    dual_obj = - kappa.dotProduct( kappa ) * 0.25  - kappa.dotProduct( b_s) - mu.dotProduct( u_s);
    
    gap = r.dotProduct( r ) - dual_obj;
}

void real_log ( DenseVector<ValueType> result, const DenseVector<ValueType> x )
{
           result[i] = log( x[i] )
        else
           result[i] = -inf;
    }
}

int main( int, char** )
{

    // number of rays

    const IndexType nRay = 50;

    // yMax x zMax is size of the 2D domain

    const IndexType yMax = 42;
    const IndexType zMax = 60;

    scai::hmemo::ContextPtr context  = scai::hmemo::Context::getContextPtr();

    // matrix D: containing the lengths of each ray in each cell
    //           (using sparse format as one ray crosses only some cells

    scai::lama::CSRSparseMatrix<ValueType> D( "input/matrix_D.mtx" );

    std::cout << "D = " << D << std::endl;

    SCAI_ASSERT_EQ_ERROR( D.getNumRows(), nRay, "size mismatch of D = " << D )
    SCAI_ASSERT_EQ_ERROR( D.getNumColumns(), yMax * zMax, "size mismatch of D = " << D )

    // T-> b

    scai::lama::DenseVector<ValueType> b( "input/vector_T.mtx" );

    std::cout << "b = " << b << std::endl;

    SCAI_ASSERT_EQ_ERROR( b.size(), nRay, "size mismatch" )

    // So -> x0, initial solution

    scai::lama::DenseVector<ValueType> x0( "input/vector_So.mtx" );

    SCAI_ASSERT_EQ_ERROR( x0.size(), yMax * zMax, "size mismatch of inital solution x0" )

    // hrz->hrz

    scai::lama::DenseVector<ValueType> hrz( "input/vector_hrz.mtx" );

    SCAI_ASSERT_EQ_ERROR( hrz.size(), mMax );

    //setup the boundary vectors

    ValueType variation = 2.0;

    scai::lama::DenseVector<ValueType> alpha( x0 );
    scai::lama::DenseVector<ValueType> beta( x0 );

    alpha *= ( 1. - variation / 100. ); //for testing: 0.01;
    beta  *= ( 1. + variation / 100. ); //for testing: 100.0;

    // construct the smoothing matrix

    CSRSparseMatrix<ValueType> L;

    setupSmoothMatrix( L, yMax, zMax );

    A.vcat( D, L );    // A = [ D; L ]

    // Extend the T vector: T_ext = [ T; zeros( size( L, 1 ) ]

    zeroExtend( T_ext, T, L.getnumRows() );

    SCAI_ASSERT_EQ_ERROR( T_ext.size(), A.getNumRows(), "serious mismatch solution vector" )

    lsqBox( sol, A, T_ext, alpha, beta );

    DenseVector res( A * sol - T );

    std::cout << "Norm = " << res.l2Norm() << std::endl;
}
