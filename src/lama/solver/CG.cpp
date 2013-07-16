/**
 * @file CG.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#include <lama/solver/CG.hpp>

// others
#include <lama/DenseVector.hpp>
#include <lama/tracing.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CG::logger, "Solver.IterativeSolver.CG" )

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
    LAMA_REGION( "Solver.CG.initialize" )
    IterativeSolver::initialize( coefficients );
    CGRuntime& runtime = getRuntime();

    runtime.mPScalar = 0.0;

    switch ( coefficients.getValueType() )
    {
    case Scalar::FLOAT:
    {
        runtime.mP.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        runtime.mQ.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        runtime.mZ.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        break;
    }
    case Scalar::DOUBLE:
    {
        runtime.mP.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        runtime.mQ.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        runtime.mZ.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        break;
    }
    default:
    {
        LAMA_THROWEXCEPTION( "Unsupported ValueType " << coefficients.getValueType() )
    }
    }

    // 'force' vector operations to be computed at the same location where coefficients reside
    runtime.mP->setContext( coefficients.getContextPtr() );
    runtime.mQ->setContext( coefficients.getContextPtr() );
    runtime.mZ->setContext( coefficients.getContextPtr() );
}

void CG::iterate()
{
    LAMA_REGION( "Solver.CG.iterate" )
    CGRuntime& runtime = getRuntime();
    Scalar lastPScalar( runtime.mPScalar );
    Scalar& pScalar = runtime.mPScalar;
    Scalar alpha;
    Scalar beta;

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
    LAMA_LOG_INFO( logger, "Doing preconditioning." )

    //CG implementation start
    if ( !mPreconditioner )
    {
        z = residual;
    }
    else
    {
        z = 0.0;
        mPreconditioner->solve( z, residual );
    }

    LAMA_LOG_INFO( logger, "Calculating pScalar." )
    pScalar = residual.dotProduct( z );
    LAMA_LOG_DEBUG( logger, "pScalar = " << pScalar )
    LAMA_LOG_INFO( logger, "Calculating p." )

    if ( this->getIterationCount() == 0 )
    {
        p = z;
    }
    else
    {
        if ( lastPScalar.getValue<double>() == 0.0 )
        {
            beta = 0.0;
        }
        else
        {
            beta = pScalar / lastPScalar;
        }

        LAMA_LOG_DEBUG( logger, "beta = " << beta )
        p = z + beta * p;
    }

    {
        LAMA_REGION( "Solver.CG.calc_q" )
        LAMA_LOG_INFO( logger, "Calculating q." )
        q = A * p;
    }

    LAMA_LOG_INFO( logger, "Calculating pqProd." )
    const Scalar pqProd = p.dotProduct( q );
    LAMA_LOG_DEBUG( logger, "pqProd = " << pqProd )

    if ( pqProd.getValue<double>() == 0.0 )
    {
        alpha = 0.0;
    }
    else
    {
        alpha = pScalar / pqProd;
    }

    LAMA_LOG_DEBUG( logger, "alpha = " << alpha )
    {
        LAMA_LOG_INFO( logger, "Calculating x." )
        LAMA_REGION( "Solver.CG.update_x" )
        x = x + alpha * p;
    }
    {
        LAMA_LOG_INFO( logger, "Updating residual." )
        LAMA_REGION( "Solver.CG.update_res" )
        residual = residual - alpha * q;
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

} //namespace lama
