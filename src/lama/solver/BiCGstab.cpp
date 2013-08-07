/**
 * @file BiCGstab.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCGstab.cpp
 * @author lschubert
 * @date 06.08.2013
 * $Id$
 */

// hpp
#include <lama/solver/BiCGstab.hpp>

// others
#include <lama/DenseVector.hpp>
#include <lama/tracing.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/Scalar.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( BiCGstab::logger, "Solver.IterativeSolver.BiCGstab" )

BiCGstab::BiCGstab( const std::string& id )
: CG( id )
{
}

BiCGstab::BiCGstab( const std::string& id, LoggerPtr logger )
: CG( id, logger )
{
}

BiCGstab::BiCGstab( const BiCGstab& other )
: CG( other )
{
}

BiCGstab::~BiCGstab()
{
}

BiCGstab::BiCGstabRuntime::BiCGstabRuntime()
    : CGRuntime()
{
}

BiCGstab::BiCGstabRuntime::~BiCGstabRuntime()
{
}

void BiCGstab::initialize( const Matrix& coefficients )
{
    LAMA_REGION( "Solver.BiCGstab.initialize" )
    CG::initialize( coefficients );
    BiCGstabRuntime& runtime = getRuntime();

    runtime.mOmega = 1.0;

    switch ( coefficients.getValueType() )
    {
    case Scalar::FLOAT:
    {
        runtime.mRes0.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        runtime.mS.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        runtime.mT.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        runtime.mTmp.reset( new DenseVector<float>( coefficients.getDistributionPtr() ) );
        break;
    }
    case Scalar::DOUBLE:
    {
        runtime.mRes0.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        runtime.mS.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        runtime.mT.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        runtime.mTmp.reset( new DenseVector<double>( coefficients.getDistributionPtr() ) );
        break;
    }
    default:
    {
        LAMA_THROWEXCEPTION( "Unsupported ValueType " << coefficients.getValueType() )
    }
    }

    // 'force' vector operations to be computed at the same location where coefficients reside
    runtime.mRes0->setContext( coefficients.getContextPtr() );
    runtime.mS->setContext( coefficients.getContextPtr() );
    runtime.mT->setContext( coefficients.getContextPtr() );
    runtime.mTmp->setContext( coefficients.getContextPtr() );
}

void BiCGstab::iterate()
{
    LAMA_REGION( "Solver.BiCGstab.iterate" )

    BiCGstabRuntime& runtime = getRuntime();
    Scalar lastPScalar( runtime.mPScalar );
    Scalar& pScalar = runtime.mPScalar;
    Scalar& lastOmega( runtime.mOmega );
    Scalar& omega = runtime.mOmega;
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
    Vector& res0 = *runtime.mRes0;

    Vector& s = *runtime.mS;
    Vector& t = *runtime.mT;
    Vector& tmp = *runtime.mTmp;

    LAMA_LOG_INFO( logger, "Doing preconditioning." )

    //BiCGstab implementation start
    if ( !mPreconditioner )
    {
        z = residual;
    }
    else
    {
//        z = 0.0;
//        mPreconditioner->solve( z, residual );
        //TODO
        LAMA_THROWEXCEPTION( "Preconditioning BiCGstab not implemented yet." )
    }


    if ( this->getIterationCount() == 0 )
    {
        pScalar = alpha = omega; // = 1
        p = res0 = z;
    }
        LAMA_LOG_INFO( logger, "Calculating pScalar." )
        pScalar = res0.dotProduct( z );
        LAMA_LOG_DEBUG( logger, "pScalar = " << pScalar )

        LAMA_LOG_INFO( logger, "Calculating p." )

        if ( lastPScalar.getValue<double>() == 0.0 )
        {
            beta = 0.0;
        }
        else
        {
            beta = ( pScalar / lastPScalar ) * ( alpha / lastOmega );
        }

        LAMA_LOG_DEBUG( logger, "beta = " << beta )
        tmp = p - lastOmega * q;
        p = z + beta * tmp;
        LAMA_LOG_TRACE( logger, "l2Norm( p ) = " << p.l2Norm() )

    {
        LAMA_REGION( "Solver.BiCGstab.calc_q" )
        LAMA_LOG_INFO( logger, "Calculating q." )
        q = A * p;
        LAMA_LOG_TRACE( logger, "l2Norm( q ) = " << q.l2Norm() )
    }

    LAMA_LOG_INFO( logger, "Calculating pqProd." )
    const Scalar pqProd = res0.dotProduct( q );
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
        LAMA_LOG_DEBUG( logger, "Calculating s.")
        LAMA_REGION( "Solver.BiCGstab.update_s")
        s = z - alpha * q;
        LAMA_LOG_TRACE( logger, "l2Norm( s ) = " << s.l2Norm() )
    }

    {
        LAMA_LOG_DEBUG( logger, "Calculating t.")
        LAMA_REGION( "Solver.BiCGstab.update_t")
        t = A * s;
        LAMA_LOG_TRACE( logger, "l2Norm( t ) = " << t.l2Norm() )
    }

    {
        LAMA_REGION( "Solver.BiCGstab.update_omega" )
        omega = t.dotProduct(s) / t.dotProduct(t);
        LAMA_LOG_TRACE( logger, "omega = " << omega )
    }

    {
        LAMA_LOG_INFO( logger, "Calculating x." )
        LAMA_REGION( "Solver.BiCGstab.update_x" )
        tmp = omega * s;
        tmp += alpha * p;
        x = x + tmp;
        LAMA_LOG_TRACE( logger, "l2Norm( x ) = " << x.l2Norm() )
    }
    {
        LAMA_LOG_INFO( logger, "Updating residual." )
        LAMA_REGION( "Solver.BiCGstab.update_res" )
        residual = s - omega * t;
        LAMA_LOG_TRACE( logger, "l2Norm( residual ) = " << residual.l2Norm() )
    }
    //BiCGstab implementation end

    mBiCGstabRuntime.mSolution.setDirty( false );
}

SolverPtr BiCGstab::copy()
{
    return SolverPtr( new BiCGstab( *this ) );
}

BiCGstab::BiCGstabRuntime& BiCGstab::getRuntime()
{
    return mBiCGstabRuntime;
}

const BiCGstab::BiCGstabRuntime& BiCGstab::getConstRuntime() const
{
    return mBiCGstabRuntime;
}

} // namespace lama
