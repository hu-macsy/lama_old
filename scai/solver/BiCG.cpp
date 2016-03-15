/**
 * @file BiCG.cpp
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
 * @brief BiCG.cpp
 * @author Lauretta Schubert
 * @date 03.07.2013
 * @since 1.1.0
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

    // runtime.mEps = std::numeric_limits<double>::epsilon() * 3;                  //CAREFUL: No abstract type
    
    switch(coefficients.getValueType()){
        case common::scalar::FLOAT:
            runtime.mEps = std::numeric_limits<float>::epsilon()*3;
            break;
        case common::scalar::DOUBLE:
            runtime.mEps = std::numeric_limits<double>::epsilon()*3;
            break;
        case common::scalar::LONG_DOUBLE:
            runtime.mEps = std::numeric_limits<long double>::epsilon()*3;
            break;
        case common::scalar::COMPLEX:
            runtime.mEps = std::numeric_limits<float>::epsilon()*3;
            break;
        case common::scalar::DOUBLE_COMPLEX:
            runtime.mEps = std::numeric_limits<double>::epsilon()*3;
            break;
        case common::scalar::LONG_DOUBLE_COMPLEX:
            runtime.mEps = std::numeric_limits<long double>::epsilon()*3;
            break;    
        default:
        SCAI_LOG_INFO(logger,"Valuetype not supported");
        break;
    }



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
    Scalar alpha;
    if( this->getIterationCount() == 0 )
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
    Scalar& eps = runtime.mEps;

    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    //BiCG implementation start
    if( !mPreconditioner )
    {
        z = residual;
        z2 = residual2;
    }
    else
    {
        z = Scalar(0.0);
        mPreconditioner->solve( z, residual );
        z2 = Scalar(0.0);
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
        Scalar beta = Scalar( 0.0 );


        if ( abs(lastPScalar) > eps )
        {
            beta = pScalar / lastPScalar;
        }

        SCAI_LOG_DEBUG( logger, "beta = " << beta << ", conj( beta ) = " << conj( beta ) )
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

    if ( abs(pqProd) < eps )
    {
         alpha = Scalar( 0.0 );
    }
    else
    {
        alpha = pScalar / pqProd;
    }

    SCAI_LOG_DEBUG( logger, "alpha = " << alpha )
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

void BiCG::print(lama::Vector& vec, size_t n){
    std::cout<<"\n";
    for(size_t i=0;i<n;++i)
        std::cout<<vec(i)<<" ";
    std::cout<<"\n";
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
