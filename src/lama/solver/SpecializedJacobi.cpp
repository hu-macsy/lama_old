/**
 * @file SpecializedJacobi.cpp
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
 * @brief CPUCSRJacobi.cpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

// hpp
#include <lama/solver/SpecializedJacobi.hpp>

// others
#include <lama/NoSyncToken.hpp>
#include <lama/Context.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/LAMAArrayUtils.hpp>

// tracing
#include <lama/tracing.hpp>

#include <boost/bind.hpp>

namespace lama
{

SpecializedJacobi::SpecializedJacobi( const std::string& id )
    : OmegaSolver( id )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, logger )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string & id, Scalar omega )
    : OmegaSolver( id, omega )
{
}

SpecializedJacobi::SpecializedJacobi( const std::string& id, Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

SpecializedJacobi::SpecializedJacobi( const SpecializedJacobi& other )
    : OmegaSolver( other )
{
}

SpecializedJacobi::SpecializedJacobiRuntime::SpecializedJacobiRuntime()
    : OmegaSolverRuntime(), mOldSolution( NULL ), mDiagonal( NULL )
{
}

SpecializedJacobi::~SpecializedJacobi()
{
    LAMA_LOG_INFO( logger, "~SpecializedJacobi" )
}

SpecializedJacobi::SpecializedJacobiRuntime::~SpecializedJacobiRuntime()
{
    LAMA_LOG_INFO( logger, "~SpecializedJacobiRuntime" )
}

void SpecializedJacobi::initialize( const Matrix& coefficients )
{
    {
        const _SparseMatrix* sparseMatrix = dynamic_cast<const _SparseMatrix*>( &coefficients );

        if ( !sparseMatrix )
        {
            LAMA_THROWEXCEPTION(
                "Coefficients matrix " << typeid(coefficients).name() << "(" << coefficients << ") is of unsupported type for SpecializedJacobi specialization (must be SparseMatrix)." );
        }
    }

    OmegaSolver::initialize( coefficients );

    LAMA_LOG_DEBUG( logger, "Initialization started" )

    LAMA_LOG_DEBUG( logger, "diagonal property of coefficients: " << coefficients.hasDiagonalProperty() )

    LAMA_ASSERT( coefficients.hasDiagonalProperty(), "Diagonal Property not set." )

    SpecializedJacobiRuntime& runtime = getRuntime();

    if ( !runtime.mOldSolution.get() )
    {
        LAMA_LOG_DEBUG( logger, "Creating old solution vector using properties of the coefficient matrix. " )
        runtime.mOldSolution.reset( Vector::createVector( runtime.mCoefficients->getValueType(),
                                    runtime.mCoefficients->getDistributionPtr() ) );
    }

    const SparseMatrix<double>* sparseDoubleCoefficients =
        dynamic_cast<const SparseMatrix<double>*>( runtime.mCoefficients );
    if ( sparseDoubleCoefficients )
    {
        LAMA_LOG_DEBUG( logger, "Creating double precison diagonal. " )
        if ( !runtime.mDiagonal.get() )
            runtime.mDiagonal.reset( _LAMAArray::create( Scalar::DOUBLE) );
        sparseDoubleCoefficients->getLocalStorage().getDiagonal( *runtime.mDiagonal );
        return;
    }
    const SparseMatrix<float>* sparseFloatCoefficients =
      dynamic_cast<const SparseMatrix<float>*>( runtime.mCoefficients );
    if ( sparseFloatCoefficients )
    {
        LAMA_LOG_DEBUG( logger, "Creating single precison diagonal. " )
        if ( !runtime.mDiagonal.get() )
            runtime.mDiagonal.reset( _LAMAArray::create( Scalar::FLOAT) );
        sparseFloatCoefficients->getLocalStorage().getDiagonal( *runtime.mDiagonal );
        return;
    }

    // has already been check in initialize, but in any case

    LAMA_THROWEXCEPTION(
        getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<T> supported)." )

//    mPointerOldSolution = &mOldSolution; --> in every solve-call
}

void SpecializedJacobi::solve( Vector& solution, const Vector& rhs )
{
    if ( getConstRuntime().mSolveInit )
    {
        LAMA_LOG_DEBUG( logger, "Previous initialization of solver found! Will be overriden!" )
    }
    solveInit( solution, rhs );
    solveImpl();
    solveFinalize();
}

void SpecializedJacobi::solveInit( Vector& solution, const Vector& rhs )
{
    //Check if oldSolution already exists, if not create copy of solution
    if ( !getConstRuntime().mOldSolution.get() )
    {
        getRuntime().mOldSolution.reset( solution.create() );

        if ( getConstRuntime().mCoefficients->getNumColumns() != getConstRuntime().mOldSolution->size() )
        {
            LAMA_THROWEXCEPTION(
                "Size of old solution vector " << *getConstRuntime().mOldSolution << " does not match number of columns of the coefficient matrix " << getConstRuntime().mCoefficients->getNumColumns() );
        }

        if ( getConstRuntime().mCoefficients->getColDistribution()
                != getConstRuntime().mOldSolution->getDistribution() )
        {
            LAMA_THROWEXCEPTION(
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
    if ( getConstRuntime().mIterations % 2 )
    {
        LAMA_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *getRuntime().mProxyOldSolution = *getRuntime().mSolution;
    }
    LAMA_LOG_DEBUG( logger, " end solve " )
}

void SpecializedJacobi::iterate()
{
    LAMA_REGION( "Solver.SpJacobi.iterate" )

    const SparseMatrix<double>* sparseDoubleCoefficients =
        dynamic_cast<const SparseMatrix<double>*>( getRuntime().mCoefficients );
    if ( sparseDoubleCoefficients )
    {
        iterateTyped( *sparseDoubleCoefficients );
        return;
    }
    const SparseMatrix<float>* sparseFloatCoefficients =
        dynamic_cast<const SparseMatrix<float>*>( getRuntime().mCoefficients );
    if ( sparseFloatCoefficients )
    {
        iterateTyped( *sparseFloatCoefficients );
        return;
    }

    // has already been check in initialize, but in any case

    LAMA_THROWEXCEPTION(
        getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<T> supported)." )
}

template<typename ValueType>
void SpecializedJacobi::iterateTyped( const SparseMatrix<ValueType>& coefficients )
{
    LAMA_REGION( "Solver.SpJacobi.iterateTyped" )

    LAMA_LOG_INFO( logger,
                   *getConstRuntime().mSolution << " = " << coefficients << " * " << *getConstRuntime().mOldSolution << " = " << *getConstRuntime().mRhs )

    if ( coefficients.getNumRows() == 0 )
    {
        LAMA_LOG_WARN( logger, "Zero sized matrix given. Won't execute any calculations in this iteration. " )
        return;
    }

    LAMA_LOG_INFO( logger, "Swap old solution and solution pointer." )

    Vector* ptr_OldSolution = &( *getRuntime().mProxyOldSolution );
    Vector* ptr_solution = &( *getRuntime().mSolution );

    getRuntime().mProxyOldSolution = ptr_solution;
    getRuntime().mSolution = ptr_OldSolution;

    //swap end now m_proxOldSolution holds the solution of the last iteration
    //and m_solution will be the output of the current iteration

    const Vector& oldSolution = getRuntime().mProxyOldSolution.getConstReference();

    //1. Check if all Vectors are DenseVectors
    if ( typeid(DenseVector<ValueType> ) == typeid( oldSolution )
            && typeid( *getRuntime().mSolution ) == typeid( oldSolution )
            && typeid( *getRuntime().mRhs ) == typeid( oldSolution ) )
    {
        LAMA_LOG_INFO( logger, "All types have the same value type." )
        const DenseVector<ValueType>& denseOldSolution = dynamic_cast<const DenseVector<ValueType>&>( oldSolution );
        DenseVector<ValueType>& denseSolution = dynamic_cast<DenseVector<ValueType>&>( *getRuntime().mSolution );
        const DenseVector<ValueType>& denseRhs = dynamic_cast<const DenseVector<ValueType>&>( *getRuntime().mRhs );

        ContextPtr localContext = coefficients.getLocalStorage().getContextPtr();

        if ( mContext )
        {
            localContext = mContext;
        }

        const ValueType omega = mOmega.getValue<ValueType>();

        // from rhs and solution we need only the local parts as LAMA arrays

        const LAMAArray<ValueType>& localRhs  = denseRhs.getLocalValues();
        LAMAArray<ValueType>& localSolution   = denseSolution.getLocalValues();
        const LAMAArray<ValueType>& localOldSolution = denseOldSolution.getLocalValues();
        LAMAArray<ValueType>&       haloOldSolution  = denseOldSolution.getHaloValues();

        const LAMAArray<ValueType>* diagonal = dynamic_cast<const LAMAArray<ValueType>*>( getRuntime().mDiagonal.get() );
 
        void ( lama::MatrixStorage<ValueType>::*jacobiIterateHalo ) ( 
            LAMAArray<ValueType>& localSolution,
            const LAMAArray<ValueType>* localDiagonal,
            const LAMAArray<ValueType>& oldHaloSolution,
            const ValueType omega ) const = &MatrixStorage<ValueType>::jacobiStepHalo;

        // will call jacobiIterateHalo( haloMatrix, localSolution, diagonal, haloOldSolution, omega )

        boost::function <void( const MatrixStorage<ValueType>* haloMatrix,
                               LAMAArray<ValueType>& localResult,
                               const LAMAArray<ValueType>& haloX )> haloF =

           boost::bind( jacobiIterateHalo, _1, _2, diagonal, _3, omega );

        if ( Matrix::SYNCHRONOUS == coefficients.getCommunicationKind() )
        {
            // For the local operation a jacobi step is done

            void ( lama::MatrixStorage<ValueType>::*jacobiStep ) ( 
                LAMAArray<ValueType>& solution,
                const LAMAArray<ValueType>& oldSolution,
                const LAMAArray<ValueType>& rhs,
                const ValueType omega ) const = &MatrixStorage<ValueType>::jacobiStep;

            // Bind the additional arguments like localRhs and omega

            boost::function <void( const MatrixStorage<ValueType>* haloMatrix,
                                   LAMAArray<ValueType>& localResult,
                                   const LAMAArray<ValueType>& localX )> localF =

                boost::bind( jacobiStep, _1, _2, _3, boost::cref( localRhs ), omega );

            coefficients.haloOperationSync( localSolution, localOldSolution, haloOldSolution, localF, haloF );
        }
        else
        {
            SyncToken* ( lama::MatrixStorage<ValueType>::*jacobiStepAsync ) ( 
                LAMAArray<ValueType>& solution,
                const LAMAArray<ValueType>& oldSolution,
                const LAMAArray<ValueType>& rhs,
                const ValueType omega ) const = &MatrixStorage<ValueType>::jacobiStepAsync;

            boost::function <SyncToken*( const MatrixStorage<ValueType>* haloMatrix,
                                         LAMAArray<ValueType>& localResult,
                                         const LAMAArray<ValueType>& localX )> localAsyncF =

                boost::bind( jacobiStepAsync, _1, _2, _3, boost::cref( localRhs ), omega );

            coefficients.haloOperationAsync( localSolution, localOldSolution, haloOldSolution, localAsyncF, haloF );
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Different types of required vectors." )
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

} // namespace lama
