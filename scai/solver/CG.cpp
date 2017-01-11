/**
 * @file CG.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief CG.cpp
 * @author Jiri Kraus
 * @date 24.08.2011
 */

// hpp
#include <scai/solver/CG.hpp>

// internal scai libraries
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( CG::logger, "Solver.IterativeSolver.CG" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

CG::CG( const std::string& id )
    : IterativeSolver( id )
{
}

CG::CG( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger )
{
}

CG::CG( const CG& other )
    : IterativeSolver( other )
{
}

CG::~CG()
{
}

CG::CGRuntime::CGRuntime()
    : IterativeSolverRuntime(), mPScalar( 0.0 )
{
}

CG::CGRuntime::~CGRuntime()
{
}

void CG::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.CG.initialize" )
    IterativeSolver::initialize( coefficients );
    CGRuntime& runtime = getRuntime();
    runtime.mPScalar = 0.0;
    runtime.mP.reset( coefficients.newDenseVector() );
    runtime.mQ.reset( coefficients.newDenseVector() );
    runtime.mZ.reset( coefficients.newDenseVector() );
}

void CG::iterate()
{
    SCAI_REGION( "Solver.CG.iterate" )
    CGRuntime& runtime = getRuntime();
    Scalar lastPScalar( runtime.mPScalar );
    Scalar& pScalar = runtime.mPScalar;

    if ( this->getIterationCount() == 0 )
    {
        this->getResidual();
    }

    Vector& residual = *runtime.mResidual;
    const Matrix& A = *runtime.mCoefficients;
    Vector& x = *runtime.mSolution;
    Vector& p = *runtime.mP;
    Vector& q = *runtime.mQ;
    Vector& z = *runtime.mZ;
    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    //CG implementation start
    if ( !mPreconditioner )
    {
        SCAI_REGION( "Solver.CG.setZ" )
        z = residual;
    }
    else
    {
        SCAI_REGION( "Solver.CG.solvePreconditioner" )
        z = Scalar( 0.0 );
        mPreconditioner->solve( z, residual );
    }

    SCAI_LOG_INFO( logger, "Calculating pScalar." )
    pScalar = z.dotProduct( residual );
    SCAI_LOG_DEBUG( logger, "pScalar = " << pScalar )
    SCAI_LOG_INFO( logger, "Calculating p." )

    if ( this->getIterationCount() == 0 )
    {
        p = z;
    }
    else
    {
        SCAI_REGION( "Solver.CG.setP" )

        // Note: lastPScalar can be very close to 0, e.g. 1e-100, is okay if pScalar is 1e-98

        Scalar beta = pScalar / lastPScalar;

        if ( Scalar( 0 ) == beta )
        {
            // ToDo: solver should terminate

            SCAI_LOG_INFO( logger, "beta = 0, can stop" )

            pScalar = lastPScalar;  // restore old value,otherwise division by zero in next step

            return;
        }

        SCAI_LOG_DEBUG( logger, "beta = " << beta << ", is p = " << pScalar << " / p_old = " << lastPScalar )

        p = z + beta * p;

        SCAI_LOG_TRACE( logger, "l2Norm( p ) = " << p.l2Norm() )
    }

    {
        SCAI_REGION( "Solver.CG.calc_q" )
        SCAI_LOG_INFO( logger, "Calculating q." )
        q = A * p;
        SCAI_LOG_TRACE( logger, "l2Norm( q ) = " << q.l2Norm() )
    }

    SCAI_LOG_INFO( logger, "Calculating pqProd." )
    const Scalar pqProd = q.dotProduct( p );
    SCAI_LOG_DEBUG( logger, "pqProd = " << pqProd )

/*    if ( pqProd == Scalar( 0.0 ) )
    {
        COMMON_THROWEXCEPTION( "Diverging due to indefinite matrix. You might try another start solution, better an adequate solver." )
    }*/

    Scalar alpha = pScalar / pqProd;

    SCAI_LOG_DEBUG( logger, "alpha = " << alpha << ", is p = " << pScalar << " / pq = " << pqProd )

    {
        SCAI_LOG_INFO( logger, "Calculating x." )
        SCAI_REGION( "Solver.CG.update_x" )
        x = x + alpha * p;
        SCAI_LOG_TRACE( logger, "l2Norm( x ) = " << x.l2Norm() )
    }
    {
        SCAI_LOG_INFO( logger, "Updating residual." )
        SCAI_REGION( "Solver.CG.update_res" )
        residual = residual - alpha * q;
        SCAI_LOG_TRACE( logger, "l2Norm( residual ) = " << residual.l2Norm() )
    }
    //CG implementation end
    mCGRuntime.mSolution.setDirty( false );
}

SolverPtr CG::copy()
{
    return SolverPtr( new CG( *this ) );
}

CG::CGRuntime& CG::getRuntime()
{
    return mCGRuntime;
}

const CG::CGRuntime& CG::getConstRuntime() const
{
    return mCGRuntime;
}

std::string CG::createValue()
{
    return "CG";
}

Solver* CG::create( const std::string name )
{
    return new CG( name );
}

void CG::writeAt( std::ostream& stream ) const
{
    stream << "CG ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
