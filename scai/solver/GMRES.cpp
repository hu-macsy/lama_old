/**
 * @file GMRES.cpp
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
 * @brief GMRES.cpp
 * @author Malte Förster
 * @date 10.04.2012
 */

// hpp
#include <scai/solver/GMRES.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/Vector.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// tracing
#include <scai/tracing.hpp>

// common
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

using utilskernel::LAMAKernel;

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, GMRES<ValueType>::logger, "Solver.IterativeSolver.GMRES" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* GMRES<ValueType>::create()
{
    return new GMRES<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType GMRES<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "GMRES" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
GMRES<ValueType>::GMRES( const std::string& id ) : 

    IterativeSolver<ValueType>( id ), 
    mKrylovDim( 10 )
{
}

template<typename ValueType>
GMRES<ValueType>::GMRES( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id, logger ), 
    mKrylovDim( 10 )
{
}

template<typename ValueType>
GMRES<ValueType>::GMRES( const GMRES<ValueType>& other ) : 

    IterativeSolver<ValueType>( other ), 
    mKrylovDim( other.mKrylovDim )
{
}

template<typename ValueType>
GMRES<ValueType>::~GMRES()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void GMRES<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.GMRES.initialize" )

    IterativeSolver<ValueType>::initialize( coefficients );

    GMRESRuntime& runtime = getRuntime();

    runtime.mCC.resize( mKrylovDim + 1 );
    runtime.mSS.resize( mKrylovDim + 1 );
    runtime.mG.resize( mKrylovDim + 1 );
    runtime.mY.resize( mKrylovDim + 1 );
    runtime.mH.resize( mKrylovDim * ( mKrylovDim + 1 ) / 2 );
    runtime.mHd.resize( mKrylovDim );

    runtime.mV.resize( mKrylovDim + 1 );   // allocate with NULL pointerso
 
    for ( IndexType i = 0; i <= mKrylovDim; ++i )
    {
        runtime.mV[i].reset( coefficients.newTargetVector() );
    }

    // 'force' vector operations to be computed at the same location where coefficients reside

    runtime.mW.reset( coefficients.newTargetVector() );
    runtime.mT.reset( coefficients.newTargetVector() );
    runtime.mX0.reset( coefficients.newTargetVector() );
}

template<typename ValueType>
void GMRES<ValueType>::setKrylovDim( IndexType krylovDim )
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

template<typename ValueType>
void GMRES<ValueType>::iterate()
{
    SCAI_REGION( "Solver.GMRES.iterate" )

    GMRESRuntime& runtime = getRuntime();

    IndexType krylovIndex = this->getIterationCount() % mKrylovDim;
    IndexType hIdxStart = krylovIndex * ( krylovIndex + 1 ) / 2;
    IndexType hIdxDiag = hIdxStart + krylovIndex;

    SCAI_LOG_INFO( logger, "GMRES( krylov dim = " << mKrylovDim << " ): iter = " << this->getIterationCount()
                           << ", inner step = " << krylovIndex )

    Vector<ValueType>& vCurrent = *runtime.mV[krylovIndex];

    const Matrix<ValueType>& A = *runtime.mCoefficients;

    // initialize in case of GMRES start/restart

    if ( krylovIndex == 0 )
    {
        SCAI_REGION( "Solver.GMRES.restartInit" )

        // Compute r0=b-Ax0

        this->getResidual();

        Vector<ValueType>& residual = *runtime.mResidual;

        // store old solution

        *runtime.mX0 = runtime.mSolution.getConstReference();

        // set first search direction vCurrent

        SCAI_LOG_DEBUG( logger, "Doing initial preconditioning." )

        if ( !mPreconditioner )
        {
            SCAI_REGION( "Solver.GMRES.setVCurrent" )
            vCurrent = residual;
        }
        else
        {
            SCAI_REGION( "Solver.GMRES.start.solvePreconditioner" )

            vCurrent.setSameValue( residual.getDistributionPtr(), 0 );

            mPreconditioner->solve( vCurrent, residual );
        }

        // normalize vCurrent

        runtime.mG[0] = vCurrent.l2Norm();

        ValueType scal = ValueType( 1 ) / runtime.mG[0];

        SCAI_LOG_DEBUG( logger, "Normalizing vCurrent with start residual " << runtime.mG[0] << "." )

        vCurrent = scal * vCurrent;
    }

    // precondition next search direction

    Vector<ValueType>& w   = *runtime.mW;
    Vector<ValueType>& tmp = *runtime.mT;

    if ( !mPreconditioner )
    {
        w = A * vCurrent;
    }
    else
    {
        SCAI_REGION( "Solver.GMRES.solvePreconditioner" )

        SCAI_LOG_DEBUG( logger, "Apply preconditioner." )

        tmp = A * vCurrent;

        w = 0;    // Note: w has been allocated with the right size

        mPreconditioner->solve( w, tmp );
    }

    // orthogonalization loop

    SCAI_LOG_DEBUG( logger, "Orthogonalization of vCurrent." )

    for ( IndexType k = 0; k <= krylovIndex; ++k )
    {
        SCAI_REGION( "Solver.GMRES.orthogonalization" )

        const Vector<ValueType>& Vk = *runtime.mV[k];

        runtime.mH[hIdxStart + k] = w.dotProduct( Vk );
        w = w - runtime.mH[hIdxStart + k] * Vk;
    }

    runtime.mHd[krylovIndex] = w.l2Norm();

    // normalize/store w in vNext (not needed in last step? Storage?)

    SCAI_LOG_DEBUG( logger, "Normalizing vNext." )

    Vector<ValueType>& vNext = *runtime.mV[krylovIndex + 1];

    ValueType scal = ValueType( 1 ) / runtime.mHd[krylovIndex];

    vNext = scal * w;

    // apply Givens rotations to new column

    SCAI_LOG_DEBUG( logger, "Apply Givens rotations." )

    for ( IndexType k = 0; k < krylovIndex; ++k )
    {
        SCAI_REGION( "Solver.GMRES.applyRotations" )

        ValueType tmp1 = runtime.mH[hIdxStart + k];
        ValueType tmp2 = runtime.mH[hIdxStart + k + 1];

        //for complex valuetype:
 
        ValueType tmpV =  common::Math::conj( runtime.mCC[k] ) * tmp1 + 
                          common::Math::conj( runtime.mSS[k] ) * tmp2;

        runtime.mH[hIdxStart + k] = tmpV;
        runtime.mH[hIdxStart + k + 1] = runtime.mCC[k] * tmp2 - runtime.mSS[k] * tmp1;
    }

    // compute new rotation

    {
        SCAI_REGION( "Solver.GMRES.computeNextRotation" )
        SCAI_LOG_DEBUG( logger, "Compute next plane rotation." )

        ValueType tmp = common::Math::sqrt(    runtime.mH[hIdxDiag] * runtime.mH[hIdxDiag]
                                            + runtime.mHd[krylovIndex] * runtime.mHd[krylovIndex] );

        runtime.mCC[krylovIndex] = runtime.mH[hIdxDiag] / tmp;
        runtime.mSS[krylovIndex] = runtime.mHd[krylovIndex] / tmp;
    }

    // update Hessenberg-system
    {
        SCAI_REGION( "Solver.GMRES.updateHessenbergSystem" )
        SCAI_LOG_DEBUG( logger, "Update Hessenberg-System." )
        runtime.mG[krylovIndex + 1] = ValueType( -1 ) * runtime.mSS[krylovIndex] * runtime.mG[krylovIndex];
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
}

template<typename ValueType>
void GMRES<ValueType>::updateX( IndexType i )
{
    SCAI_REGION( "Solver.GMRES.updateX" )

    // back-substitution Hessenberg system H * y = g
    // H stored in column 'packed' order

    SCAI_LOG_INFO( logger, "updating solution within krylov dimensions i+1 = " << i + 1 )

    GMRESRuntime& runtime = getRuntime();

    // implementation using LAPACK

    for ( IndexType j = 0; j <= i; ++j )
    {
        runtime.mY[j] = runtime.mG[j];
    }

    // ContextPtr context = getCoefficients().getContextPtr();


    hmemo::ContextPtr context = hmemo::Context::getHostPtr();
    static LAMAKernel<blaskernel::BLASKernelTrait::tptrs<ValueType> > tptrs;

    const IndexType n = i + 1;   // size of matrix 
    const IndexType nrhs = 1;

    tptrs[context]( CblasUpper, common::MatrixOp::NORMAL, CblasNonUnit, n, nrhs, &runtime.mH[0],
                    &runtime.mY[0], n );

    // Update of solution vector

    Vector<ValueType>& x = runtime.mSolution.getReference();  // -> dirty

    // reset x to x0

    if ( i != 0 )
    {
        x = *runtime.mX0;
    }

    // update x
    // TODO: Add linar combination method

    for ( IndexType k = 0; k <= i; ++k )
    {
        x += runtime.mY[k] * ( *runtime.mV[k]);   // axpy calls
    }
}


/* ========================================================================= */
/*       Getter runtime                                                      */
/* ========================================================================= */

template<typename ValueType>
typename GMRES<ValueType>::GMRESRuntime& GMRES<ValueType>::getRuntime()
{
    return mGMRESRuntime;
}

template<typename ValueType>
const typename GMRES<ValueType>::GMRESRuntime& GMRES<ValueType>::getRuntime() const
{
    return mGMRESRuntime;
}

/* ========================================================================= */
/*       virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
GMRES<ValueType>* GMRES<ValueType>::copy()
{
    return new GMRES<ValueType>( *this );
}

template<typename ValueType>
void GMRES<ValueType>::writeAt( std::ostream& stream ) const
{
    const char* typeId = common::TypeTraits<ValueType>::id();

    stream << "GMRES<" << typeId << "> ( id = " << this->getId() << ", krylov dim = " << mKrylovDim
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( GMRES, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
