/**
 * @file Jacobi.cpp
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
 * @brief Implementation of class Jacobi.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/Jacobi.hpp>

// local library
#include <scai/utilskernel/HArrayUtils.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/macros/loop.hpp>

namespace scai
{

using tasking::SyncToken;

using utilskernel::HArrayUtils;
using lama::Matrix;
using lama::Vector;
using lama::DenseVector;
using lama::Scalar;

namespace solver
{

Jacobi::Jacobi( const std::string& id )
    : OmegaSolver( id )
{
}

Jacobi::Jacobi( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, logger )
{
}

Jacobi::Jacobi( const std::string& id, lama::Scalar omega )
    : OmegaSolver( id, omega )
{
}

Jacobi::Jacobi( const std::string& id, lama::Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

Jacobi::Jacobi( const Jacobi& other )
    : OmegaSolver( other )
{
}

Jacobi::JacobiRuntime::JacobiRuntime()
    : OmegaSolverRuntime(), mOldSolution(), mDiagonal()
{
}

Jacobi::~Jacobi()
{
    SCAI_LOG_INFO( logger, "~Jacobi" )
}

Jacobi::JacobiRuntime::~JacobiRuntime()
{
    SCAI_LOG_INFO( logger, "~JacobiRuntime" )
}

void Jacobi::initialize( const Matrix& coefficients )
{
    using hmemo::_HArray;

    if ( coefficients.getMatrixKind() == Matrix::DENSE )
    {
        COMMON_THROWEXCEPTION(
            "Coefficients matrix " << typeid( coefficients ).name() << "(" << coefficients << ") is of unsupported type for Jacobi specialization (must be SparseMatrix)." );
    }

    OmegaSolver::initialize( coefficients );
    SCAI_LOG_DEBUG( logger, "Initialization started" )
    SCAI_LOG_DEBUG( logger, "diagonal property of coefficients: " << coefficients.hasDiagonalProperty() )
    SCAI_ASSERT( coefficients.hasDiagonalProperty(), "Diagonal Property not set." )
    JacobiRuntime& runtime = getRuntime();

    if ( !runtime.mOldSolution.get() )
    {
        SCAI_LOG_DEBUG( logger, "Creating old solution vector using properties of the coefficient matrix. " )
        runtime.mOldSolution.reset( runtime.mCoefficients->newDenseVector() );
    }

    if ( runtime.mCoefficients->getMatrixKind() == Matrix::SPARSE )
    {
        if ( !runtime.mDiagonal.get() )
        {
            runtime.mDiagonal.reset( _HArray::create( runtime.mCoefficients->getValueType() ) );
        }

        runtime.mCoefficients->getLocalStorage().getDiagonal( *runtime.mDiagonal );
    }
    else
    {
        COMMON_THROWEXCEPTION    (
            getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<ValueType> supported)." )
    }

//    mPointerOldSolution = &mOldSolution; --> in every solve-call
}

void Jacobi::solveInit( Vector& solution, const Vector& rhs )
{
    //Check if oldSolution already exists, if not create copy of solution
    if ( !getConstRuntime().mOldSolution.get() )
    {
        // Important: method newVector creats vector with same context as solution

        getRuntime().mOldSolution.reset( solution.newVector() );

        if ( getConstRuntime().mCoefficients->getNumColumns() != getConstRuntime().mOldSolution->size() )
        {
            COMMON_THROWEXCEPTION(
                "Size of old solution vector " << *getConstRuntime().mOldSolution << " does not match number of columns of the coefficient matrix " << getConstRuntime().mCoefficients->getNumColumns() );
        }

        if ( getConstRuntime().mCoefficients->getColDistribution() != getConstRuntime().mOldSolution->getDistribution() )
        {
            COMMON_THROWEXCEPTION(
                "Distribution of " << *getConstRuntime().mOldSolution << " = " << getConstRuntime().mOldSolution->getDistribution() << " does not match column distribution of " << *getConstRuntime().mCoefficients << " = " << getConstRuntime().mCoefficients->getColDistribution() );
        }
    }

    getRuntime().mProxyOldSolution = getConstRuntime().mOldSolution.get();
    IterativeSolver::solveInit( solution, rhs );
}

void Jacobi::solveFinalize()
{
//    MF: ?????
//    if( &( mProxyOldSolution.getConstReference() ) ==
//        &( mSolution.getConstReference() ) )
    if ( getConstRuntime().mIterations % 2 )
    {
        SCAI_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *getRuntime().mProxyOldSolution = *getRuntime().mSolution;
    }

    SCAI_LOG_DEBUG( logger, " end solve " )
}

void Jacobi::iterate()
{
    SCAI_REGION( "Solver.SpJacobi.iterate" )
    // for each supported arithmetic type we have to dynamic cast and instatiate typed version
#define SCAI_SOLVER_TYPE_CAST( _type )                                                                                  \
    {                                                                                               \
        const lama::SparseMatrix<_type>* sparseTypedCoefficients =                                  \
                dynamic_cast<const lama::SparseMatrix<_type>*>( getRuntime().mCoefficients );       \
        if ( sparseTypedCoefficients )                                                              \
        {                                                                                           \
            iterateTyped( *sparseTypedCoefficients );                                               \
            return;                                                                                 \
        }                                                                                           \
    }
    SCAI_COMMON_LOOP( SCAI_SOLVER_TYPE_CAST, SCAI_NUMERIC_TYPES_HOST )
#undef SCAI_SOLVER_TYPE_CAST
    // has already been check in initialize, but in any case
    COMMON_THROWEXCEPTION        (
        getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<ValueType> supported)." )
}

template<typename ValueType>
void Jacobi::iterateTyped( const lama::SparseMatrix<ValueType>& coefficients )
{
    using hmemo::HArray;
    SCAI_REGION( "Solver.SpJacobi.iterateTyped" )
    SCAI_LOG_INFO( logger,
                   *getConstRuntime().mSolution << " = " << coefficients << " * " << *getConstRuntime().mOldSolution << " = " << *getConstRuntime().mRhs )

    if ( coefficients.getNumRows() == 0 )
    {
        SCAI_LOG_WARN( logger, "Zero sized matrix given. Won't execute any calculations in this iteration. " )
        return;
    }

    SCAI_LOG_INFO( logger, "Swap old solution and solution pointer." )
    Vector* ptr_OldSolution = &( *getRuntime().mProxyOldSolution );
    Vector* ptr_solution = &( *getRuntime().mSolution );
    getRuntime().mProxyOldSolution = ptr_solution;
    getRuntime().mSolution = ptr_OldSolution;
    //swap end now m_proxOldSolution holds the solution of the last iteration
    //and m_solution will be the output of the current iteration
    const Vector& oldSolution = getRuntime().mProxyOldSolution.getConstReference();

    //1. Check if all Vectors are DenseVectors
    if (  DenseVector<ValueType>::createValue() == oldSolution.getCreateValue()
            && ( *getRuntime().mSolution ).getCreateValue() == oldSolution.getCreateValue()
            && ( *getRuntime().mRhs ).getCreateValue() == oldSolution.getCreateValue() )
//    if( typeid(DenseVector<ValueType> ) == typeid( oldSolution )
//            && typeid( *getRuntime().mSolution ) == typeid( oldSolution )
//            && typeid( *getRuntime().mRhs ) == typeid( oldSolution ) )
    {
        SCAI_LOG_INFO( logger, "All types have the same value type." )
        const DenseVector<ValueType>& denseOldSolution = dynamic_cast<const DenseVector<ValueType>&>( oldSolution );
        DenseVector<ValueType>& denseSolution = dynamic_cast<DenseVector<ValueType>&>( *getRuntime().mSolution );
        const DenseVector<ValueType>& denseRhs = dynamic_cast<const DenseVector<ValueType>&>( *getRuntime().mRhs );
        hmemo::ContextPtr localContext = coefficients.getLocalStorage().getContextPtr();

        if ( mContext )
        {
            localContext = mContext;
        }

        const ValueType omega = mOmega.getValue<ValueType>();

        // from rhs and solution we need only the local parts as LAMA arrays
        const HArray<ValueType>& localRhs = denseRhs.getLocalValues();

        HArray<ValueType>& localSolution = denseSolution.getLocalValues();

        const HArray<ValueType>& localOldSolution = denseOldSolution.getLocalValues();

        HArray<ValueType>& haloOldSolution = denseOldSolution.getHaloValues();

        const HArray<ValueType>* diagonal = dynamic_cast<const HArray<ValueType>*>( getRuntime().mDiagonal.get() );

        using namespace scai::common;  // placeholders are also needed

        void ( scai::lama::MatrixStorage<ValueType>::*jacobiIterateHalo )(
            HArray<ValueType>& localSolution,
            const HArray<ValueType>& localDiagonal,
            const HArray<ValueType>& oldHaloSolution,
            const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateHalo;

        // will call jacobiIterateHalo( haloMatrix, localSolution, diagonal, haloOldSolution, omega )
        function <
        void(
            const lama::MatrixStorage<ValueType>* haloMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& haloX ) > haloF =
                bind( jacobiIterateHalo, _1, _2, cref( *diagonal ), _3, omega );

        if ( Matrix::SYNCHRONOUS == coefficients.getCommunicationKind() )
        {
            // For the local operation a jacobi step is done
            void ( lama::MatrixStorage<ValueType>::*jacobiIterate )(
                HArray<ValueType>& solution,
                const HArray<ValueType>& oldSolution,
                const HArray<ValueType>& rhs,
                const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterate;
            // Bind the additional arguments like localRhs and omega
            function <
            void(
                const lama::MatrixStorage<ValueType>* haloMatrix,
                HArray<ValueType>& localResult,
                const HArray<ValueType>& localX ) > localF =
                    bind( jacobiIterate, _1, _2, _3, cref( localRhs ), omega );
            coefficients.haloOperationSync( localSolution, localOldSolution, haloOldSolution, localF, haloF );
        }
        else
        {
            SyncToken* ( lama::MatrixStorage<ValueType>::*jacobiIterateAsync )(
                HArray<ValueType>& solution,
                const HArray<ValueType>& oldSolution,
                const HArray<ValueType>& rhs,
                const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateAsync;
            function <
            SyncToken*(
                const lama::MatrixStorage<ValueType>* haloMatrix,
                HArray<ValueType>& localResult,
                const HArray<ValueType>& localX ) > localAsyncF =
                    bind( jacobiIterateAsync, _1, _2, _3, cref( localRhs ), omega );
            coefficients.haloOperationAsync( localSolution, localOldSolution, haloOldSolution, localAsyncF, haloF );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "Different types of required vectors." )
    }
}

Jacobi::JacobiRuntime& Jacobi::getRuntime()
{
    return mJacobiRuntime;
}

const Jacobi::JacobiRuntime& Jacobi::getConstRuntime() const
{
    return mJacobiRuntime;
}

SolverPtr Jacobi::copy()
{
    return SolverPtr( new Jacobi( *this ) );
}

void Jacobi::writeAt( std::ostream& stream ) const
{
    stream << "Jacobi ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string Jacobi::createValue()
{
    return "Jacobi";
}

Solver* Jacobi::create( const std::string name )
{
    return new Jacobi( name );
}

} /* end namespace solver */

} /* end namespace scai */
