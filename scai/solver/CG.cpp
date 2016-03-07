/**
 * @file CG.cpp
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
 * @brief CG.cpp
 * @author Jiri Kraus
 * @date 24.08.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/CG.hpp>

// internal scai libraries
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>

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
    
    runtime.mEps = std::numeric_limits<double>::epsilon() * 3;                  //CAREFUL: No abstract type

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
