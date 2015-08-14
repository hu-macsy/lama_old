/**
 * @file GMRES.cpp
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
 * @brief GMRES.cpp
 * @author Malte Förster
 * @date 10.04.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/GMRES.hpp>

// others
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/LAMAInterface.hpp>

// tracing
#include <scai/tracing.hpp>

#include <omp.h>

using common::unique_ptr;
using common::scoped_array;

namespace lama
{

LAMA_LOG_DEF_LOGGER( GMRES::logger, "Solver.IterativeSolver.GMRES" )

GMRES::GMRES( const std::string& id )
    : IterativeSolver( id ), mKrylovDim( 10 ), totalIterationTime( 0.0 ), totalPreconditionerTime( 0.0 )
{
}

GMRES::GMRES( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger ), mKrylovDim( 10 )
{
}

GMRES::GMRES( const GMRES& other )
    : IterativeSolver( other ), mKrylovDim( other.mKrylovDim )
{
}

GMRES::GMRESRuntime::GMRESRuntime()
    : IterativeSolverRuntime(), mCC(new double[1]), mSS(new double[1]), mG(new double[1]), mY(new double[1]), mH(new double[1]), mHd(new double[1]), mV( 0 ), mW( 0 ),
        mT( 0 ), mX0( 0 )
{
}

GMRES::~GMRES()
{
}

double GMRES::getAverageIterationTime() const
{
    return ( ( this->totalIterationTime - this->totalPreconditionerTime ) / this->getIterationCount() );
}

double GMRES::getAveragePreconditionerTime() const
{
    return ( this->totalPreconditionerTime / this->getIterationCount() );
}
GMRES::GMRESRuntime::~GMRESRuntime()
{
    if( mV != 0 )
    {
        for( unsigned int i = 0; i < mV->size(); ++i )
        {
            delete ( *mV )[i];
        }

        delete mV;
    }

    mV = 0;

    if( mW != 0 )
    {
        delete mW;
    }

    mW = 0;

    if( mT != 0 )
    {
        delete mT;
    }

    mT = 0;

    if( mX0 != 0 )
    {
        delete mX0;
    }

    mX0 = 0;
}

void GMRES::initialize( const Matrix& coefficients )
{
    LAMA_REGION( "Solver.GMRES.initialize" )

    IterativeSolver::initialize( coefficients );

    GMRESRuntime& runtime = getRuntime();

    if( runtime.mV != 0 )
    {
        for( unsigned int i = 0; i < runtime.mV->size(); ++i )
        {
            delete ( *runtime.mV )[i];
        }

        delete runtime.mV;
    }

    runtime.mV = 0;

    if( runtime.mW != 0 )
    {
        delete runtime.mW;
    }

    runtime.mW = 0;

    if( runtime.mT != 0 )
    {
        delete runtime.mT;
    }

    runtime.mT = 0;

    if( runtime.mX0 != 0 )
    {
        delete runtime.mX0;
    }

    runtime.mX0 = 0;

    scoped_array<double>& mCC = runtime.mCC;
    scoped_array<double>& mSS = runtime.mSS;
    scoped_array<double>& mG = runtime.mG;
    scoped_array<double>& mY = runtime.mY;
    scoped_array<double>& mH = runtime.mH;
    scoped_array<double>& mHd = runtime.mHd;

    mCC.reset( new double[mKrylovDim + 1] );
    mSS.reset( new double[mKrylovDim + 1] );
    mG.reset( new double[mKrylovDim + 1] );
    mY.reset( new double[mKrylovDim + 1] );

    mH.reset( new double[( mKrylovDim * ( mKrylovDim + 1 ) ) / 2] );
    mHd.reset( new double[mKrylovDim] );

    runtime.mV = new std::vector<Vector*>( mKrylovDim + 1, 0 );

    switch( coefficients.getValueType() )
    {
        case common::scalar::FLOAT:
        {

            ( *runtime.mV )[0] = new DenseVector<float>( coefficients.getDistributionPtr() );

            runtime.mW = new DenseVector<float>( coefficients.getDistributionPtr() );
            runtime.mT = new DenseVector<float>( coefficients.getDistributionPtr() );
            runtime.mX0 = new DenseVector<float>( coefficients.getDistributionPtr() );

            break;
        }

        case common::scalar::DOUBLE:
        {
            ( *runtime.mV )[0] = new DenseVector<double>( coefficients.getDistributionPtr() );

            runtime.mW = new DenseVector<double>( coefficients.getDistributionPtr() );
            runtime.mT = new DenseVector<double>( coefficients.getDistributionPtr() );
            runtime.mX0 = new DenseVector<double>( coefficients.getDistributionPtr() );

            break;
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Unsupported ValueType " << coefficients.getValueType() )
        }
    }

    // 'force' vector operations to be computed at the same location where coefficients reside

    ( *runtime.mV )[0]->setContext( coefficients.getContextPtr() );
    runtime.mW->setContext( coefficients.getContextPtr() );
    runtime.mT->setContext( coefficients.getContextPtr() );
    runtime.mX0->setContext( coefficients.getContextPtr() );

    totalIterationTime = 0.0;
    totalPreconditionerTime = 0.0;
}

void GMRES::setKrylovDim( unsigned int krylovDim )
{
    LAMA_LOG_DEBUG( logger, " Krylov dimension set to " << krylovDim )
    mKrylovDim = krylovDim;
}

/*
 * 1. Compute r0=b-Ax0
 * 1*. Compute r0=b-Ax0, prec.solve(v1,r0)
 * 1**. Compute r0=b-Ax0
 * normalize to v1
 *
 * Define H(m+1,m)=0
 *
 * 2. loop j=1,m
 * 3.   wj=Avj
 * 3*.  prec.solve(wj,Avj)
 * 3**. prec.solve(tmp,vj), wj=A*tmp
 * 4.   loop i=1,j
 * 5.       scalar product
 * 6.       orthogonalize
 * 7.   normalize
 * 8.   if (norm==0) break loop
 * 9. minimize ||beta*e1-Hy||
 * 10.xm=x0+Vm*ym
 * 10*.xm=x0+Vm*ym
 * 10**. prec.solve(tmp, Vm*ym), xm=x0+tmp
 */

void GMRES::iterate()
{
    LAMA_REGION( "Solver.GMRES.iterate" )

    GMRESRuntime& runtime = getRuntime();

    double iterationTimeStart = omp_get_wtime();

    unsigned int krylovIndex = this->getIterationCount() % mKrylovDim;
    unsigned int hIdxStart = krylovIndex * ( krylovIndex + 1 ) / 2;
    unsigned int hIdxDiag = hIdxStart + krylovIndex;
    LAMA_LOG_INFO( logger, "GMRES("<<mKrylovDim<<"): Inner Step "<<krylovIndex<<"." )
    Vector& vCurrent = *( ( *runtime.mV )[krylovIndex] );
    const Matrix& A = ( *runtime.mCoefficients );

    // lazy allocation structure mV
    if( !( *runtime.mV )[krylovIndex + 1] )
    {
        LAMA_REGION( "Solver.GMRES.setMV" )

        switch( A.getValueType() )
        {
            case common::scalar::FLOAT:
            {
                ( *runtime.mV )[krylovIndex + 1] = new DenseVector<float>( A.getDistributionPtr() );
                break;
            }

            case common::scalar::DOUBLE:
            {
                ( *runtime.mV )[krylovIndex + 1] = new DenseVector<double>( A.getDistributionPtr() );
                break;
            }

            default:
            {
                COMMON_THROWEXCEPTION( "Unsupported ValueType " << A.getValueType() )
            }
        }

        ( *runtime.mV )[krylovIndex + 1]->setContext( A.getContextPtr() );
    }

    // initialize in case of GMRES start/restart
    if( krylovIndex == 0 )
    {
        LAMA_REGION( "Solver.GMRES.restartInit" )
        // Compute r0=b-Ax0
        this->getResidual();
        Vector& residual = ( *runtime.mResidual );

        // store old solution
        *runtime.mX0 = runtime.mSolution.getConstReference();

        // set first search direction vCurrent
        LAMA_LOG_INFO( logger, "Doing initial preconditioning." )

        if( !mPreconditioner )
        {
            LAMA_REGION( "Solver.GMRES.setVCurrent" )
            vCurrent = residual;
        }
        else
        {
            LAMA_REGION( "Solver.GMRES.start.solvePreconditioner" )
            vCurrent = 0.0;
            double preconditionerTimeStart = omp_get_wtime();
            mPreconditioner->solve( vCurrent, residual );
            totalPreconditionerTime += omp_get_wtime() - preconditionerTimeStart;
        }

        // normalize vCurrent
        runtime.mG[0] = vCurrent.l2Norm().getValue<double>();
        double scal = 1.0 / runtime.mG[0];
        LAMA_LOG_DEBUG( logger, "Normalizing vCurrent with start residual "<< runtime.mG[0] <<"." )
        vCurrent = scal * vCurrent;
    }

    // precondition next search direction
    Vector& w = ( *runtime.mW );
    Vector& tmp = ( *runtime.mT );

    LAMA_LOG_INFO( logger, "Doing preconditioning." )

    if( !mPreconditioner )
    {
        w = A * vCurrent;
    }
    else
    {
        LAMA_REGION( "Solver.GMRES.solvePreconditioner" )
        tmp = A * vCurrent;
        w = 0.0;
        double preconditionerTimeStart = omp_get_wtime();
        mPreconditioner->solve( w, tmp );
        totalPreconditionerTime += omp_get_wtime() - preconditionerTimeStart;
    }

    // orthogonalization loop
    LAMA_LOG_DEBUG( logger, "Orthogonalization of vCurrent." )

    for( unsigned int k = 0; k <= krylovIndex; ++k )
    {
        LAMA_REGION( "Solver.GMRES.orthogonalization" )
        const Vector& Vk = *( ( *runtime.mV )[k] );
        runtime.mH[hIdxStart + k] = ( w.dotProduct( Vk ) ).getValue<double>();
        w = w - runtime.mH[hIdxStart + k] * Vk;
    }

    runtime.mHd[krylovIndex] = w.l2Norm().getValue<double>();

    // normalize/store w in vNext (not needed in last step? Storage?)
    LAMA_LOG_DEBUG( logger, "Normalizing vNext." )
    Vector& vNext = *( *runtime.mV )[krylovIndex + 1];
    double scal = 1.0 / runtime.mHd[krylovIndex];
    vNext = scal * w;

    // apply Givens rotations to new column
    LAMA_LOG_DEBUG( logger, "Apply Givens rotations." )

    for( unsigned int k = 0; k < krylovIndex; ++k )
    {
        LAMA_REGION( "Solver.GMRES.applyRotations" )
        double tmp1 = runtime.mH[hIdxStart + k];
        double tmp2 = runtime.mH[hIdxStart + k + 1];
        runtime.mH[hIdxStart + k] = runtime.mCC[k] * tmp1 + runtime.mSS[k] * tmp2;
        runtime.mH[hIdxStart + k + 1] = runtime.mCC[k] * tmp2 - runtime.mSS[k] * tmp1;
    }

    // compute new rotation
    {
        LAMA_REGION( "Solver.GMRES.computeNextRotation" )
        LAMA_LOG_DEBUG( logger, "Compute next plane rotation." )
        double tmp = std::sqrt(
                         runtime.mH[hIdxDiag] * runtime.mH[hIdxDiag]
                         + runtime.mHd[krylovIndex] * runtime.mHd[krylovIndex] );
        runtime.mCC[krylovIndex] = runtime.mH[hIdxDiag] / tmp;
        runtime.mSS[krylovIndex] = runtime.mHd[krylovIndex] / tmp;
    }

    // update Hessenberg-system
    {
        LAMA_REGION( "Solver.GMRES.updateHessenbergSystem" )
        LAMA_LOG_DEBUG( logger, "Update Hessenberg-System." )
        runtime.mG[krylovIndex + 1] = -1.0 * runtime.mSS[krylovIndex] * runtime.mG[krylovIndex];
        runtime.mG[krylovIndex] = runtime.mCC[krylovIndex] * runtime.mG[krylovIndex];
        runtime.mH[hIdxDiag] = runtime.mCC[krylovIndex] * runtime.mH[hIdxDiag]
                               + runtime.mSS[krylovIndex] * runtime.mHd[krylovIndex];
        LAMA_LOG_DEBUG( logger, "New Residual estimate "<< abs(runtime.mG[krylovIndex+1]) << "." )
    }

    // do (partial) update to solution (currently very expensive to do every iteration)
    // TODO do this more efficiently? Only do update if residual evaluation
    // required or krylov-subspace completely filled
    //if (krylovIndex == mKrylovDim-1)
    updateX( krylovIndex );

    totalIterationTime += omp_get_wtime() - iterationTimeStart;
}

void GMRES::updateX( unsigned int i )
{
    LAMA_REGION( "Solver.GMRES.updateX" )

    // back-substitution Hessenberg system H*y=g
    // H stored in column 'packed' order
    LAMA_LOG_DEBUG( logger, "Updating X within krylov dimensions i+1 = " << i+1 )

    GMRESRuntime& runtime = getRuntime();

    // implementation using LAPACK
    for( unsigned int j = 0; j <= i; ++j )
    {
        runtime.mY[j] = runtime.mG[j];
    }

    // ContextPtr context = getCoefficients().getContextPtr();

    memory::ContextPtr context = memory::Context::getContextPtr( memory::context::Host );

    LAMA_INTERFACE_FN_t( tptrs, context, BLAS, LAPACK, double );

    int info = tptrs( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, i + 1, 1, runtime.mH.get(),
                      runtime.mY.get(), i + 1 );

    LAMA_LOG_DEBUG( logger, "tptrs returned with code = " << info )

    // Update of solution vector
    Vector& x = runtime.mSolution.getReference();

    // reset x to x0
    if( i != 0 )
    {
        x = *runtime.mX0;
    }

    // update x
    // TODO: Add linar combination method
    for( unsigned int k = 0; k <= i; ++k )
    {
        const Vector& Vk = *( ( *runtime.mV )[k] );
        x = x + runtime.mY[k] * Vk;
    }
}

GMRES::GMRESRuntime& GMRES::getRuntime()
{
    return mGMRESRuntime;
}

const GMRES::GMRESRuntime& GMRES::getConstRuntime() const
{
    return mGMRESRuntime;
}

SolverPtr GMRES::copy()
{
    return SolverPtr( new GMRES( *this ) );
}

} //namespace lama
