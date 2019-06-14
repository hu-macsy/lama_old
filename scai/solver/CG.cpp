/**
 * @file CG.cpp
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
 * @brief CG.cpp
 * @author Jiri Kraus
 * @date 24.08.2011
 */

// hpp
#include <scai/solver/CG.hpp>

// internal scai libraries
#include <scai/lama/Vector.hpp>
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, CG<ValueType>::logger, "Solver.IterativeSolver.CG" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* CG<ValueType>::create()
{
    return new CG<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType CG<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "CG" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
CG<ValueType>::CG( const std::string& id ) : 

    IterativeSolver<ValueType>( id )
{
}

template<typename ValueType>
CG<ValueType>::CG( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
CG<ValueType>::CG( const CG& other ) : 

    IterativeSolver<ValueType>( other )
{
}

template<typename ValueType>
CG<ValueType>::~CG()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void CG<ValueType>::CGRuntime::initialize()
{
    SCAI_ASSERT_ERROR( mCoefficients, "initialize of CGRuntime, but no coefficients set" )

    mPScalar = ValueType( 0 );

    // allocate runtime vectors p, q, z with target space + context of matrix

    mP.reset( mCoefficients->newTargetVector() );
    mQ.reset( mCoefficients->newTargetVector() );
    mZ.reset( mCoefficients->newTargetVector() );
}

template<typename ValueType>
void CG<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.CG.initialize" )

    IterativeSolver<ValueType>::initialize( coefficients );

    getRuntime().initialize();
}

/* ========================================================================= */
/*    solve : one iteration                                                  */
/* ========================================================================= */

template<typename ValueType>
void CG<ValueType>::iterate()
{
    SCAI_LOG_INFO( logger, "CG.iterate, iter = " << IterativeSolver<ValueType>::getIterationCount() )

    SCAI_REGION( "Solver.CG.iterate" )
    CGRuntime& runtime = getRuntime();

    ValueType& pScalar    = runtime.mPScalar;   // keep ref for vallue p that is kept between it
    ValueType lastPScalar = pScalar;            // copy of last p 

    if ( IterativeSolver<ValueType>::getIterationCount() == 0 )
    {
        this->getResidual();
    }

    Vector<ValueType>& residual = *runtime.mResidual;

    const Matrix<ValueType>& A = *runtime.mCoefficients;  // coefficient matrix is pointer

    Vector<ValueType>& x = runtime.mSolution.getReference();  // will be updated
    Vector<ValueType>& p = *runtime.mP;
    Vector<ValueType>& q = *runtime.mQ;
    Vector<ValueType>& z = *runtime.mZ;

    SCAI_LOG_INFO( logger, "Doing preconditioning." )

    // CG implementation start
    if ( !this->mPreconditioner )
    {
        SCAI_REGION( "Solver.CG.setZ" )
        z = residual;
    }
    else
    {
        SCAI_REGION( "Solver.CG.solvePreconditioner" )
        z.setSameValue( A.getRowDistributionPtr(), 0 );
        this->mPreconditioner->solve( z, residual );
    }

    SCAI_LOG_INFO( logger, "Calculating pScalar." )
    pScalar = z.dotProduct( residual );
    SCAI_LOG_DEBUG( logger, "pScalar = " << pScalar )
    SCAI_LOG_INFO( logger, "Calculating p." )

    if ( IterativeSolver<ValueType>::getIterationCount() == 0 )
    {
        p = z;
    }
    else
    {
        SCAI_REGION( "Solver.CG.setP" )

        // Note: lastPScalar can be very close to 0, e.g. 1e-100, is okay if pScalar is 1e-98

        ValueType beta = pScalar / lastPScalar;

        if ( ValueType( 0 ) == beta )
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
    const ValueType pqProd = q.dotProduct( p );
    SCAI_LOG_DEBUG( logger, "pqProd = " << pqProd )

    /*    if ( pqProd == Scalar( 0.0 ) )
    {
        COMMON_THROWEXCEPTION( "Diverging due to indefinite matrix. You might try another start solution, better an adequate solver." )
    }*/

    ValueType alpha = pScalar / pqProd;

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

template<typename ValueType>
CG<ValueType>* CG<ValueType>::copy()
{
    return new CG<ValueType>( *this );    // copy by using the copy constructor
}

template<typename ValueType>
typename CG<ValueType>::CGRuntime& CG<ValueType>::getRuntime()
{
    return mCGRuntime;
}

template<typename ValueType>
const typename CG<ValueType>::CGRuntime& CG<ValueType>::getRuntime() const
{
    return mCGRuntime;
}

template<typename ValueType>
void CG<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "CG<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId() 
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( CG, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
