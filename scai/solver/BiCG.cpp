/**
 * @file BiCG.cpp
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
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, BiCG<ValueType>::logger, "Solver.IterativeSolver.BiCG" )

using lama::Matrix;
using lama::Vector;
using lama::DenseVector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* BiCG<ValueType>::create()
{
    return new BiCG<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType BiCG<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "BiCG" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
BiCG<ValueType>::BiCG( const std::string& id ) : 

    CG<ValueType>( id )
{
}

template<typename ValueType>
BiCG<ValueType>::BiCG( const std::string& id, LoggerPtr logger ) : 

    CG<ValueType>( id, logger )
{
}

template<typename ValueType>
BiCG<ValueType>::BiCG( const BiCG<ValueType>& other ) : 

    CG<ValueType>( other )
{
}

template<typename ValueType>
BiCG<ValueType>::~BiCG()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void BiCG<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.BiCG.initialize" )

    CG<ValueType>::initialize( coefficients );

    BiCGRuntime& runtime = getRuntime();

    runtime.mPScalar2 = ValueType( 0 );

    hmemo::ContextPtr ctx = coefficients.getContextPtr();

    runtime.mConjTransposeA.reset( coefficients.newMatrix() );
    runtime.mConjTransposeA->assignTranspose( coefficients );
    runtime.mConjTransposeA->conj();

    runtime.mP2.setContextPtr( ctx );
    runtime.mQ2.setContextPtr( ctx );
    runtime.mZ2.setContextPtr( ctx );
    runtime.mResidual2.setContextPtr( ctx );

    // temporary vectors will be allocated at their use
}

template<typename ValueType>
void BiCG<ValueType>::iterate()
{
    SCAI_REGION( "Solver.BiCG.iterate" )

    BiCGRuntime& runtime = getRuntime();

    ValueType  lastPScalar = runtime.mPScalar;
    ValueType& pScalar = runtime.mPScalar;

    if ( this->getIterationCount() == 0 )
    {
        this->getResidual();
        this->getResidual2();
    }

    DenseVector<ValueType>& residual = runtime.mResidual;
    DenseVector<ValueType>& residual2 = runtime.mResidual2;

    const Matrix<ValueType>& A = *runtime.mCoefficients;
    const Matrix<ValueType>& Act = *runtime.mConjTransposeA;  

    Vector<ValueType>& x = runtime.mSolution.getReference(); // ->dirty

    DenseVector<ValueType>& p = runtime.mP;
    DenseVector<ValueType>& p2 = runtime.mP2;
    DenseVector<ValueType>& q = runtime.mQ;
    DenseVector<ValueType>& q2 = runtime.mQ2;
    DenseVector<ValueType>& z = runtime.mZ;
    DenseVector<ValueType>& z2 = runtime.mZ2;

    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    if ( !mPreconditioner )
    {
        z = residual;
        z2 = residual2;
    }
    else
    {
        z.setSameValue( residual.getDistributionPtr(), 0 );
        mPreconditioner->solve( z, residual );
        z2.setSameValue( residual2.getDistributionPtr(), 0 );
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
        ValueType beta =  pScalar / lastPScalar;

        SCAI_LOG_DEBUG( logger, "beta = " << beta << ", is p = " << pScalar << " / p_old = " << lastPScalar )

        p = z + beta * p;
        SCAI_LOG_TRACE( logger, "l2Norm( p ) = " << p.l2Norm() )
        p2 = z2 + common::Math::conj( beta ) * p2;
        SCAI_LOG_TRACE( logger, "l2Norm( p2 ) = " << p2.l2Norm() )
    }

    {
        SCAI_REGION( "Solver.BiCG.calc_q" )
        SCAI_LOG_INFO( logger, "Calculating q." )
        q = A * p;
        SCAI_LOG_TRACE( logger, "l2Norm( q ) = " << q.l2Norm() )
        q2 = Act * p2; // conjTranspose( A ) * p2 
        SCAI_LOG_TRACE( logger, "l2Norm( q2 ) = " << q2.l2Norm() )
    }

    SCAI_LOG_INFO( logger, "Calculating pqProd." )
    const ValueType pqProd = p2.dotProduct( q );
    SCAI_LOG_DEBUG( logger, "pqProd = " << pqProd )

    /*    if ( pqProd == Scalar( 0.0 ) )
    {
        COMMON_THROWEXCEPTION( "Diverging due to indefinite matrix. You might try another start solution, better an adequate solver." )
    }*/

    ValueType alpha = pScalar / pqProd;

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
        residual2 = residual2 - common::Math::conj( alpha ) * q2;
        //residual2 = residual2 - alpha * q2;
        SCAI_LOG_TRACE( logger, "l2Norm( residual2 ) = " << residual.l2Norm() )
    }
    //BiCG implementation end
    mBiCGRuntime.mSolution.setDirty( false );
}

template<typename ValueType>
const Vector<ValueType>& BiCG<ValueType>::getResidual2() const
{
    SCAI_LOG_DEBUG( logger, "getResidual2 of solver " << *this )

    const BiCGRuntime& runtime = getRuntime();

    SCAI_ASSERT_DEBUG( runtime.mSolveInit, "no residual for unintialized solver." )

    SCAI_LOG_DEBUG( logger, "calculating residual of = " << runtime.mSolution.getConstReference() )

    const Vector<ValueType>& solution = runtime.mSolution.getConstReference(); // only read
    const Vector<ValueType>& rhs = *runtime.mRhs;
    const Matrix<ValueType>& Act = *runtime.mConjTransposeA;

    runtime.mResidual2 = rhs - Act * solution;   // rhs - conjTranspose( A ) * solution

    return runtime.mResidual2;
}

template<typename ValueType>
void BiCG<ValueType>::print( lama::Vector<ValueType>& vec, size_t n )
{
    std::cout << "\n";

    for ( size_t i = 0; i < n; ++i )
    {
        ValueType val = vec[i];
        std::cout << val << " ";
    }

    std::cout << "\n";
}

/* ========================================================================= */
/*       Getter runtime                                                      */
/* ========================================================================= */

template<typename ValueType>
typename BiCG<ValueType>::BiCGRuntime& BiCG<ValueType>::getRuntime()
{
    return mBiCGRuntime;
}

template<typename ValueType>
const typename BiCG<ValueType>::BiCGRuntime& BiCG<ValueType>::getRuntime() const
{
    return mBiCGRuntime;
}

/* ========================================================================= */
/*       Virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
BiCG<ValueType>* BiCG<ValueType>::copy()
{
    return new BiCG( *this );
}

template<typename ValueType>
void BiCG<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "BiCG<" << this->getValueType() << "> ( id = " << this->getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( BiCG, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
