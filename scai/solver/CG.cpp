/**
 * @file CG.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief CG.cpp
 * @author Jiri Kraus
 * @date 24.08.2011
 */

// hpp
#include <scai/solver/CG.hpp>

// local library
#include <scai/solver/mepr/SolverEps.hpp>

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
    runtime.mEps = mepr::SolverEps<SCAI_ARITHMETIC_HOST_LIST>::get( coefficients.getValueType() ) * 3.0;
    
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
    const Scalar& eps = runtime.mEps;
    Scalar alpha;
    Scalar beta;


    if( this->getIterationCount() == 0 )
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
    if( !mPreconditioner )
    {
        SCAI_REGION( "Solver.CG.setZ" )
        z = residual;
    }
    else
    {
        SCAI_REGION( "Solver.CG.solvePreconditioner" )
        z = Scalar(0.0);
        mPreconditioner->solve( z, residual );
    }

    SCAI_LOG_INFO( logger, "Calculating pScalar." )
    pScalar = z.dotProduct( residual );
    SCAI_LOG_DEBUG( logger, "pScalar = " << pScalar )
    SCAI_LOG_INFO( logger, "Calculating p." )

    if( this->getIterationCount() == 0 )
    {
        p = z;
    }
    else
    {
        SCAI_REGION( "Solver.CG.setP" )

        if( lastPScalar.getValue<double>() < eps )  //scalar is small
        {
            beta = 0.0;
        }
        else
        {
            beta = pScalar / lastPScalar;
        }

        SCAI_LOG_DEBUG( logger, "beta = " << beta )
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

    if( pqProd.getValue<double>() < eps )   //scalar is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = pScalar / pqProd;
    }

    SCAI_LOG_DEBUG( logger, "alpha = " << alpha )
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
