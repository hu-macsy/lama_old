/**
 * @file GMRES.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief GMRES.cpp
 * @author Malte Förster
 * @date 10.04.2012
 */

// hpp
#include <scai/solver/GMRES.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// tracing
#include <scai/tracing.hpp>

// common
#include <scai/common/Walltime.hpp>

namespace scai
{

namespace solver
{

using utilskernel::LAMAKernel;

using lama::Matrix;
using lama::Vector;
using lama::DenseVector;
using lama::Scalar;

using common::unique_ptr;
using common::scoped_array;

SCAI_LOG_DEF_LOGGER( GMRES::logger, "Solver.IterativeSolver.GMRES" )

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
    : IterativeSolverRuntime(), mCC( new double[1] ), mSS( new double[1] ), mG( new double[1] ), mY( new double[1] ), mH( new double[1] ), mHd( new double[1] ), mV( 0 ), mW( 0 ),
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
    if ( mV != 0 )
    {
        for ( unsigned int i = 0; i < mV->size(); ++i )
        {
            delete ( *mV )[i];
        }

        delete mV;
    }

    mV = 0;

    if ( mW != 0 )
    {
        delete mW;
    }

    mW = 0;

    if ( mT != 0 )
    {
        delete mT;
    }

    mT = 0;

    if ( mX0 != 0 )
    {
        delete mX0;
    }

    mX0 = 0;
}

void GMRES::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.GMRES.initialize" )
    IterativeSolver::initialize( coefficients );
    GMRESRuntime& runtime = getRuntime();

    if ( runtime.mV != 0 )
    {
        for ( unsigned int i = 0; i < runtime.mV->size(); ++i )
        {
            delete ( *runtime.mV )[i];
        }

        delete runtime.mV;
    }

    runtime.mV = 0;

    if ( runtime.mW != 0 )
    {
        delete runtime.mW;
    }

    runtime.mW = 0;

    if ( runtime.mT != 0 )
    {
        delete runtime.mT;
    }

    runtime.mT = 0;

    if ( runtime.mX0 != 0 )
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
    // 'force' vector operations to be computed at the same location where coefficients reside
    ( *runtime.mV )[0] = coefficients.newDenseVector();
    runtime.mW = coefficients.newDenseVector();
    runtime.mT = coefficients.newDenseVector();
    runtime.mX0 = coefficients.newDenseVector();
    totalIterationTime = 0.0;
    totalPreconditionerTime = 0.0;
}

void GMRES::setKrylovDim( unsigned int krylovDim )
{
    SCAI_LOG_DEBUG( logger, " Krylov dimension set to " << krylovDim )
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
    SCAI_REGION( "Solver.GMRES.iterate" )
    GMRESRuntime& runtime = getRuntime();
    double iterationTimeStart = common::Walltime::get();
    unsigned int krylovIndex = this->getIterationCount() % mKrylovDim;
    unsigned int hIdxStart = krylovIndex * ( krylovIndex + 1 ) / 2;
    unsigned int hIdxDiag = hIdxStart + krylovIndex;
    SCAI_LOG_INFO( logger, "GMRES(" << mKrylovDim << "): Inner Step " << krylovIndex << "." )
    Vector& vCurrent = *( ( *runtime.mV )[krylovIndex] );
    const Matrix& A = ( *runtime.mCoefficients );

    // lazy allocation structure mV
    if ( !( *runtime.mV )[krylovIndex + 1] )
    {
        SCAI_REGION( "Solver.GMRES.setMV" )
        ( *runtime.mV )[krylovIndex + 1] = A.newDenseVector();
    }

    // initialize in case of GMRES start/restart
    if ( krylovIndex == 0 )
    {
        SCAI_REGION( "Solver.GMRES.restartInit" )
        // Compute r0=b-Ax0
        this->getResidual();
        Vector& residual = ( *runtime.mResidual );
        // store old solution
        *runtime.mX0 = runtime.mSolution.getConstReference();
        // set first search direction vCurrent
        SCAI_LOG_INFO( logger, "Doing initial preconditioning." )

        if ( !mPreconditioner )
        {
            SCAI_REGION( "Solver.GMRES.setVCurrent" )
            vCurrent = residual;
        }
        else
        {
            SCAI_REGION( "Solver.GMRES.start.solvePreconditioner" )
            vCurrent = 0.0;
            double preconditionerTimeStart = common::Walltime::get();
            mPreconditioner->solve( vCurrent, residual );
            totalPreconditionerTime += common::Walltime::get() - preconditionerTimeStart;
        }

        // normalize vCurrent
        runtime.mG[0] = vCurrent.l2Norm().getValue<double>();
        double scal = 1.0 / runtime.mG[0];
        SCAI_LOG_DEBUG( logger, "Normalizing vCurrent with start residual " << runtime.mG[0] << "." )
        vCurrent = scal * vCurrent;
    }

    // precondition next search direction
    Vector& w = ( *runtime.mW );
    Vector& tmp = ( *runtime.mT );
    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    if ( !mPreconditioner )
    {
        w = A * vCurrent;
    }
    else
    {
        SCAI_REGION( "Solver.GMRES.solvePreconditioner" )
        tmp = A * vCurrent;
        w = 0.0;
        double preconditionerTimeStart = common::Walltime::get();
        mPreconditioner->solve( w, tmp );
        totalPreconditionerTime += common::Walltime::get() - preconditionerTimeStart;
    }

    // orthogonalization loop
    SCAI_LOG_DEBUG( logger, "Orthogonalization of vCurrent." )

    for ( unsigned int k = 0; k <= krylovIndex; ++k )
    {
        SCAI_REGION( "Solver.GMRES.orthogonalization" )
        const Vector& Vk = *( ( *runtime.mV )[k] );
        runtime.mH[hIdxStart + k] = ( w.dotProduct( Vk ) ).getValue<double>();
        w = w - runtime.mH[hIdxStart + k] * Vk;
    }

    runtime.mHd[krylovIndex] = w.l2Norm().getValue<double>();
    // normalize/store w in vNext (not needed in last step? Storage?)
    SCAI_LOG_DEBUG( logger, "Normalizing vNext." )
    Vector& vNext = *( *runtime.mV )[krylovIndex + 1];
    double scal = 1.0 / runtime.mHd[krylovIndex];
    vNext = scal * w;
    // apply Givens rotations to new column
    SCAI_LOG_DEBUG( logger, "Apply Givens rotations." )

    for ( unsigned int k = 0; k < krylovIndex; ++k )
    {
        SCAI_REGION( "Solver.GMRES.applyRotations" )
        double tmp1 = runtime.mH[hIdxStart + k];
        double tmp2 = runtime.mH[hIdxStart + k + 1];
        //for complex valuetype:
        runtime.mH[hIdxStart + k] = lama::conj( runtime.mCC[k] ).getValue<ScalarRepType>() * tmp1 + lama::conj( runtime.mSS[k] ).getValue<ScalarRepType>() * tmp2;
        runtime.mH[hIdxStart + k + 1] = runtime.mCC[k] * tmp2 - runtime.mSS[k] * tmp1;
    }

    // compute new rotation
    {
        SCAI_REGION( "Solver.GMRES.computeNextRotation" )
        SCAI_LOG_DEBUG( logger, "Compute next plane rotation." )
        double tmp = std::sqrt(
                         runtime.mH[hIdxDiag] * runtime.mH[hIdxDiag]
                         + runtime.mHd[krylovIndex] * runtime.mHd[krylovIndex] );
        runtime.mCC[krylovIndex] = runtime.mH[hIdxDiag] / tmp;
        runtime.mSS[krylovIndex] = runtime.mHd[krylovIndex] / tmp;
    }
    // update Hessenberg-system
    {
        SCAI_REGION( "Solver.GMRES.updateHessenbergSystem" )
        SCAI_LOG_DEBUG( logger, "Update Hessenberg-System." )
        runtime.mG[krylovIndex + 1] = -1.0 * runtime.mSS[krylovIndex] * runtime.mG[krylovIndex];
        runtime.mG[krylovIndex] = runtime.mCC[krylovIndex] * runtime.mG[krylovIndex];
        runtime.mH[hIdxDiag] = runtime.mCC[krylovIndex] * runtime.mH[hIdxDiag]
                               + runtime.mSS[krylovIndex] * runtime.mHd[krylovIndex];
        SCAI_LOG_DEBUG( logger, "New Residual estimate " << common::Math::abs( runtime.mG[krylovIndex + 1] ) << "." )
    }
    // do (partial) update to solution (currently very expensive to do every iteration)
    // TODO do this more efficiently? Only do update if residual evaluation
    // required or krylov-subspace completely filled
    //if (krylovIndex == mKrylovDim-1)
    updateX( krylovIndex );
    totalIterationTime += common::Walltime::get() - iterationTimeStart;
}

void GMRES::updateX( unsigned int i )
{
    SCAI_REGION( "Solver.GMRES.updateX" )
    // back-substitution Hessenberg system H*y=g
    // H stored in column 'packed' order
    SCAI_LOG_DEBUG( logger, "Updating X within krylov dimensions i+1 = " << i + 1 )
    GMRESRuntime& runtime = getRuntime();

    // implementation using LAPACK
    for ( unsigned int j = 0; j <= i; ++j )
    {
        runtime.mY[j] = runtime.mG[j];
    }

    // ContextPtr context = getCoefficients().getContextPtr();
    hmemo::ContextPtr context = hmemo::Context::getHostPtr();
    static LAMAKernel<blaskernel::BLASKernelTrait::tptrs<double> > tptrs;

    tptrs[context]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, i + 1, 1, runtime.mH.get(),
                    runtime.mY.get(), i + 1 );

    // Update of solution vector

    Vector& x = runtime.mSolution.getReference();

    // reset x to x0
    if ( i != 0 )
    {
        x = *runtime.mX0;
    }

    // update x
    // TODO: Add linar combination method
    for ( unsigned int k = 0; k <= i; ++k )
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

void GMRES::writeAt( std::ostream& stream ) const
{
    stream << "GMRES ( id = " << mId << ", krylov dim = " << mKrylovDim
           << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string GMRES::createValue()
{
    return "GMRES";
}

Solver* GMRES::create( const std::string name )
{
    return new GMRES( name );
}

} /* end namespace solver */

} /* end namespace scai */
