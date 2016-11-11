/**
 * @file BiCG.cpp
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
 * @brief BiCG.cpp
 * @author Lauretta Schubert
 * @date 03.07.2013
 */

// hpp
#include <scai/solver/BiCG.hpp>

// internal scai libraries
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/storage/MatrixStorage.hpp>

#include <scai/tracing.hpp>
#include <scai/common/ScalarType.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( BiCG::logger, "Solver.IterativeSolver.BiCG" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

BiCG::BiCG( const std::string& id )
    : CG( id )
{
}

BiCG::BiCG( const std::string& id, LoggerPtr logger )
    : CG( id, logger )
{
}

BiCG::BiCG( const BiCG& other )
    : CG( other )
{
}

BiCG::~BiCG()
{
}

BiCG::BiCGRuntime::BiCGRuntime()
    : CGRuntime(), mPScalar2( 0.0 ), mResidual2()
{
}

BiCG::BiCGRuntime::~BiCGRuntime()
{
}

void BiCG::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.BiCG.initialize" )
    CG::initialize( coefficients );
    BiCGRuntime& runtime = getRuntime();
    runtime.mPScalar2 = 0.0;
    runtime.mTransposeA.reset( coefficients.newMatrix() );
    common::scalar::ScalarType type = coefficients.getValueType();
    runtime.mP2.reset( Vector::getDenseVector( type, coefficients.getRowDistributionPtr() ) );
    runtime.mQ2.reset( Vector::getDenseVector( type, coefficients.getRowDistributionPtr() ) );
    runtime.mZ2.reset( Vector::getDenseVector( type, coefficients.getRowDistributionPtr() ) );
    runtime.mResidual2.reset( Vector::getDenseVector( type, coefficients.getRowDistributionPtr() ) );
    runtime.mTransposeA->assignTranspose( coefficients );
    runtime.mTransposeA->conj();
    // 'force' vector operations to be computed at the same location where coefficients reside
    runtime.mP2->setContextPtr( coefficients.getContextPtr() );
    runtime.mQ2->setContextPtr( coefficients.getContextPtr() );
    runtime.mZ2->setContextPtr( coefficients.getContextPtr() );
    runtime.mResidual2->setContextPtr( coefficients.getContextPtr() );
}

void BiCG::iterate()
{
    SCAI_REGION( "Solver.BiCG.iterate" )
    BiCGRuntime& runtime = getRuntime();
    Scalar lastPScalar( runtime.mPScalar );
    Scalar& pScalar = runtime.mPScalar;

    if ( this->getIterationCount() == 0 )
    {
        this->getResidual();
        this->getResidual2();
    }

    Vector& residual = *runtime.mResidual;
    Vector& residual2 = *runtime.mResidual2;
    const Matrix& A = *runtime.mCoefficients;
    const Matrix& transA = *runtime.mTransposeA;
    Vector& x = *runtime.mSolution;
    Vector& p = *runtime.mP;
    Vector& p2 = *runtime.mP2;
    Vector& q = *runtime.mQ;
    Vector& q2 = *runtime.mQ2;
    Vector& z = *runtime.mZ;
    Vector& z2 = *runtime.mZ2;
    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    //BiCG implementation start
    if ( !mPreconditioner )
    {
        z = residual;
        z2 = residual2;
    }
    else
    {
        z = Scalar( 0.0 );
        mPreconditioner->solve( z, residual );
        z2 = Scalar( 0.0 );
// THIS IS WRONG!!
// Instead of solving P * z2 = residual2 we need to solve P^H * z2 = residual2
// where P is the preconditioner
        mPreconditioner->solve( z2, residual2 );
        // print(residual,4);
        // print(residual2,4);
    }

    SCAI_LOG_INFO( logger, "Calculating pScalar." )
    pScalar = z2.dotProduct( z );
    SCAI_LOG_DEBUG( logger, "pScalar = " << pScalar )
    SCAI_LOG_INFO( logger, "Calculating p." )

    if ( this->getIterationCount() == 0 )
    {
        p = z;
        p2 = z2;
    }
    else
    {
        Scalar beta =  pScalar / lastPScalar;

        SCAI_LOG_DEBUG( logger, "beta = " << beta << ", is p = " << pScalar << " / p_old = " << lastPScalar )

        p = z + beta * p;
        SCAI_LOG_TRACE( logger, "l2Norm( p ) = " << p.l2Norm() )
        p2 = z2 + conj( beta ) * p2;
        SCAI_LOG_TRACE( logger, "l2Norm( p2 ) = " << p2.l2Norm() )
    }

    {
        SCAI_REGION( "Solver.BiCG.calc_q" )
        SCAI_LOG_INFO( logger, "Calculating q." )
        q = A * p;
        SCAI_LOG_TRACE( logger, "l2Norm( q ) = " << q.l2Norm() )
        q2 = transA * p2; //p2 * A;
        SCAI_LOG_TRACE( logger, "l2Norm( q2 ) = " << q2.l2Norm() )
    }

    SCAI_LOG_INFO( logger, "Calculating pqProd." )
    const Scalar pqProd = p2.dotProduct( q );
    SCAI_LOG_DEBUG( logger, "pqProd = " << pqProd )

    if( pqProd == Scalar( 0.0 ) )
    {
        COMMON_THROWEXCEPTION( "Diverging due to indefinite matrix. You might try another start solution, better an adequate solver." )
    }

    Scalar alpha = pScalar / pqProd;

    SCAI_LOG_DEBUG( logger, "alpha = " << alpha << ", is p = " << pScalar << " / pq = " << pqProd )


    {
        SCAI_LOG_INFO( logger, "Calculating x." )
        SCAI_REGION( "Solver.BiCG.update_x" )
        x = x + alpha * p;
        SCAI_LOG_TRACE( logger, "l2Norm( x ) = " << x.l2Norm() )
    }
    {
        SCAI_LOG_INFO( logger, "Updating residual." )
        SCAI_REGION( "Solver.BiCG.update_res" )
        residual = residual - alpha * q;
        SCAI_LOG_TRACE( logger, "l2Norm( residual ) = " << residual.l2Norm() )
        residual2 = residual2 - conj( alpha ) * q2;
        //residual2 = residual2 - alpha * q2;
        SCAI_LOG_TRACE( logger, "l2Norm( residual2 ) = " << residual.l2Norm() )
    }
    //BiCG implementation end
    mBiCGRuntime.mSolution.setDirty( false );
}

const Vector& BiCG::getResidual2() const
{
    SCAI_LOG_DEBUG( logger, "getResidual2 of solver " << mId )
    const BiCGRuntime& runtime = getConstRuntime();
    SCAI_ASSERT_DEBUG( runtime.mCoefficients, "mCoefficients == NULL" )
    SCAI_ASSERT_DEBUG( runtime.mRhs, "mRhs == NULL" )
    //mLogger->logMessage(LogLevel::completeInformation,"Request for residual received.\n");
    SCAI_LOG_DEBUG( logger, "calculating residual of = " << runtime.mSolution.getConstReference() )
    //mLogger->logMessage(LogLevel::completeInformation,"Residual needs revaluation.\n");
    mLogger->startTimer( "ResidualTimer" );
    *runtime.mResidual2 = *runtime.mRhs;
    *runtime.mResidual2 -= ( *runtime.mTransposeA ) * runtime.mSolution.getConstReference() ;
    mLogger->stopTimer( "ResidualTimer" );
    mLogger->logTime( "ResidualTimer", LogLevel::completeInformation, "Revaluation of residual took [s]: " );
    mLogger->stopAndResetTimer( "ResidualTimer" );
    return ( *runtime.mResidual2 );
}

void BiCG::print( lama::Vector& vec, size_t n )
{
    std::cout << "\n";

    for ( size_t i = 0; i < n; ++i )
    {
        std::cout << vec( i ) << " ";
    }

    std::cout << "\n";
}

SolverPtr BiCG::copy()
{
    return SolverPtr( new BiCG( *this ) );
}

BiCG::BiCGRuntime& BiCG::getRuntime()
{
    return mBiCGRuntime;
}

const BiCG::BiCGRuntime& BiCG::getConstRuntime() const
{
    return mBiCGRuntime;
}

std::string BiCG::createValue()
{
    return "BiCG";
}

Solver* BiCG::create( const std::string name )
{
    return new BiCG( name );
}

void BiCG::writeAt( std::ostream& stream ) const
{
    stream << "BiCG ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
