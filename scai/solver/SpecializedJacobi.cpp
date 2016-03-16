/**
 * @file SpecializedJacobi.cpp
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
 * @brief Implementation of class SpecializedJacobi.
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/SpecializedJacobi.hpp>

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

SpecializedJacobi::SpecializedJacobi( const std::string& id )
    : OmegaSolver( id )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, logger )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string & id, lama::Scalar omega )
    : OmegaSolver( id, omega )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string& id, lama::Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

SpecializedJacobi::SpecializedJacobi( const SpecializedJacobi& other )
    : OmegaSolver( other )
{
}

SpecializedJacobi::SpecializedJacobiRuntime::SpecializedJacobiRuntime()
    : OmegaSolverRuntime(), mOldSolution(), mDiagonal()
{
}

SpecializedJacobi::~SpecializedJacobi()
{
    SCAI_LOG_INFO( logger, "~SpecializedJacobi" )
}

SpecializedJacobi::SpecializedJacobiRuntime::~SpecializedJacobiRuntime()
{
    SCAI_LOG_INFO( logger, "~SpecializedJacobiRuntime" )
}

void SpecializedJacobi::initialize( const Matrix& coefficients )
{
    using hmemo::_HArray;

    if( coefficients.getMatrixKind() == Matrix::DENSE )
    {
        COMMON_THROWEXCEPTION(
            "Coefficients matrix " << typeid(coefficients).name() << "(" << coefficients << ") is of unsupported type for SpecializedJacobi specialization (must be SparseMatrix)." );
    }

    OmegaSolver::initialize( coefficients );

    SCAI_LOG_DEBUG( logger, "Initialization started" )

    SCAI_LOG_DEBUG( logger, "diagonal property of coefficients: " << coefficients.hasDiagonalProperty() )

    SCAI_ASSERT( coefficients.hasDiagonalProperty(), "Diagonal Property not set." )

    SpecializedJacobiRuntime& runtime = getRuntime();

    if( !runtime.mOldSolution.get() )
    {
        SCAI_LOG_DEBUG( logger, "Creating old solution vector using properties of the coefficient matrix. " )
        runtime.mOldSolution.reset( runtime.mCoefficients->newDenseVector() );
    }

    if( runtime.mCoefficients->getMatrixKind() == Matrix::SPARSE )
    {
        if( !runtime.mDiagonal.get() )
        {
            runtime.mDiagonal.reset( _HArray::create( runtime.mCoefficients->getValueType()));
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

void SpecializedJacobi::solveInit( Vector& solution, const Vector& rhs )
{
    //Check if oldSolution already exists, if not create copy of solution
    if( !getConstRuntime().mOldSolution.get() )
    {
        getRuntime().mOldSolution.reset( Vector::create( solution.getCreateValue() ) );

        if( getConstRuntime().mCoefficients->getNumColumns() != getConstRuntime().mOldSolution->size() )
        {
            COMMON_THROWEXCEPTION(
                "Size of old solution vector " << *getConstRuntime().mOldSolution << " does not match number of columns of the coefficient matrix " << getConstRuntime().mCoefficients->getNumColumns() );
        }

        if( getConstRuntime().mCoefficients->getColDistribution() != getConstRuntime().mOldSolution->getDistribution() )
        {
            COMMON_THROWEXCEPTION(
                "Distribution of " << *getConstRuntime().mOldSolution << " = " << getConstRuntime().mOldSolution->getDistribution() << " does not match column distribution of " << *getConstRuntime().mCoefficients << " = " << getConstRuntime().mCoefficients->getColDistribution() );
        }
    }

    getRuntime().mProxyOldSolution = getConstRuntime().mOldSolution.get();

    IterativeSolver::solveInit( solution, rhs );
}

void SpecializedJacobi::solveFinalize()
{
//    MF: ?????
//    if( &( mProxyOldSolution.getConstReference() ) ==
//        &( mSolution.getConstReference() ) )
    if( getConstRuntime().mIterations % 2 )
    {
        SCAI_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *getRuntime().mProxyOldSolution = *getRuntime().mSolution;
    }

    SCAI_LOG_DEBUG( logger, " end solve " )
}

void SpecializedJacobi::iterate()
{
    SCAI_REGION( "Solver.SpJacobi.iterate" )

    // for each supported arithmetic type we have to dynamic cast and instatiate typed version

#define SCAI_SOLVER_TYPE_CAST( _type )                                                                                  \
    {                                                                                                             \
        const lama::SparseMatrix<_type>* sparseTypedCoefficients =                             \
                reinterpret_cast<const lama::SparseMatrix<_type>*>( getRuntime().mCoefficients );  \
        if ( sparseTypedCoefficients )                                                                            \
                   {                                                                                              \
            iterateTyped( *sparseTypedCoefficients );                                                             \
            return;                                                                                               \
        }                                                                                                         \
    }

    SCAI_COMMON_TYPELOOP( ARITHMETIC_HOST_CNT, SCAI_SOLVER_TYPE_CAST, ARITHMETIC_HOST )

#undef SCAI_SOLVER_TYPE_CAST

    // has already been check in initialize, but in any case

    COMMON_THROWEXCEPTION        (
        getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<ValueType> supported)." )
}

template<typename ValueType>
void SpecializedJacobi::iterateTyped( const lama::SparseMatrix<ValueType>& coefficients )
{
    using hmemo::HArray;

    SCAI_REGION( "Solver.SpJacobi.iterateTyped" )

    SCAI_LOG_INFO( logger,
                   *getConstRuntime().mSolution << " = " << coefficients << " * " << *getConstRuntime().mOldSolution << " = " << *getConstRuntime().mRhs )

    if( coefficients.getNumRows() == 0 )
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
    if( typeid(DenseVector<ValueType> ) == typeid( oldSolution )
            && typeid( *getRuntime().mSolution ) == typeid( oldSolution )
            && typeid( *getRuntime().mRhs ) == typeid( oldSolution ) )
    {
        SCAI_LOG_INFO( logger, "All types have the same value type." )
        const DenseVector<ValueType>& denseOldSolution = dynamic_cast<const DenseVector<ValueType>&>( oldSolution );
        DenseVector<ValueType>& denseSolution = dynamic_cast<DenseVector<ValueType>&>( *getRuntime().mSolution );
        const DenseVector<ValueType>& denseRhs = dynamic_cast<const DenseVector<ValueType>&>( *getRuntime().mRhs );

        hmemo::ContextPtr localContext = coefficients.getLocalStorage().getContextPtr();

        if( mContext )
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

        void (scai::lama::MatrixStorage<ValueType>::*jacobiIterateHalo)(
            HArray<ValueType>& localSolution,
            const HArray<ValueType>& localDiagonal,
            const HArray<ValueType>& oldHaloSolution,
            const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateHalo;

        // will call jacobiIterateHalo( haloMatrix, localSolution, diagonal, haloOldSolution, omega )

        function<
        void(
            const lama::MatrixStorage<ValueType>* haloMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& haloX )> haloF =

                bind( jacobiIterateHalo, _1, _2, cref( *diagonal ), _3, omega );

        if( Matrix::SYNCHRONOUS == coefficients.getCommunicationKind() )
        {
            // For the local operation a jacobi step is done

            void (lama::MatrixStorage<ValueType>::*jacobiIterate)(
                HArray<ValueType>& solution,
                const HArray<ValueType>& oldSolution,
                const HArray<ValueType>& rhs,
                const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterate;

            // Bind the additional arguments like localRhs and omega

            function<
            void(
                const lama::MatrixStorage<ValueType>* haloMatrix,
                HArray<ValueType>& localResult,
                const HArray<ValueType>& localX )> localF =

                    bind( jacobiIterate, _1, _2, _3, cref( localRhs ), omega );

            coefficients.haloOperationSync( localSolution, localOldSolution, haloOldSolution, localF, haloF );
        }
        else
        {
            SyncToken* (lama::MatrixStorage<ValueType>::*jacobiIterateAsync)(
                HArray<ValueType>& solution,
                const HArray<ValueType>& oldSolution,
                const HArray<ValueType>& rhs,
                const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateAsync;

            function<
            SyncToken*(
                const lama::MatrixStorage<ValueType>* haloMatrix,
                HArray<ValueType>& localResult,
                const HArray<ValueType>& localX )> localAsyncF =

                    bind( jacobiIterateAsync, _1, _2, _3, cref( localRhs ), omega );

            coefficients.haloOperationAsync( localSolution, localOldSolution, haloOldSolution, localAsyncF, haloF );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "Different types of required vectors." )
    }
}

SpecializedJacobi::SpecializedJacobiRuntime& SpecializedJacobi::getRuntime()
{
    return mSpecializedJacobiRuntime;
}

const SpecializedJacobi::SpecializedJacobiRuntime& SpecializedJacobi::getConstRuntime() const
{
    return mSpecializedJacobiRuntime;
}

SolverPtr SpecializedJacobi::copy()
{
    return SolverPtr( new SpecializedJacobi( *this ) );
}

void SpecializedJacobi::writeAt( std::ostream& stream ) const
{
    stream << "SpecializedJacobi ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string SpecializedJacobi::createValue()
{
	return "SpecializedJacobi";
}

Solver* SpecializedJacobi::create( const std::string name )
{
	return new SpecializedJacobi( name );
}

} /* end namespace solver */

} /* end namespace scai */
