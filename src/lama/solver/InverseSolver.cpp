/**
 * @file InverseSolver.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief InverseSolver.cpp
 * @author Jiri Kraus
 * @date 08.06.2011
 * $Id$
 */

// hpp
#include <lama/solver/InverseSolver.hpp>

// others
#include <lama/solver/h_frag.h>

#include <lama/norm/L2Norm.hpp>

#include <lama/matrix/DenseMatrix.hpp>

#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/tracing.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/distribution/NoDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>

// boost
#include <boost/scoped_array.hpp>

#include <sstream>

namespace lama
{

LAMA_LOG_DEF_LOGGER( InverseSolver::logger, "Solver.InverseSolver" );

InverseSolver::InverseSolver( const std::string & id )
    : Solver( id )
{
}

InverseSolver::InverseSolver( const std::string & id, LoggerPtr logger )
    : Solver( id, logger )
{
    LAMA_LOG_INFO( InverseSolver::logger, "InverseSolver, id = " << id );
}

InverseSolver::InverseSolver( const InverseSolver& other )
    : Solver( other )
{
    LAMA_LOG_INFO( InverseSolver::logger, "InverseSolver, id = " << other.mId );
}

InverseSolver::InverseSolverRuntime::InverseSolverRuntime()
    : SolverRuntime()
{
}

InverseSolver::~InverseSolver()
{
    LAMA_LOG_INFO( logger, "~InverseSolver" );
}

InverseSolver::InverseSolverRuntime::~InverseSolverRuntime()
{
}

void InverseSolver::initialize( const Matrix& coefficients )
{
    LAMA_REGION( "Solver.Inverse.intialize" );

    LAMA_LOG_INFO( logger, "Initializing with " << coefficients );

    getRuntime().mInverse = MatrixPtr( coefficients.create().release() );

    getRuntime().mInverse->invert( coefficients );

    getRuntime().mInverse->setContext( coefficients.getContextPtr() );

    getRuntime().mInverse->prefetch();

    Solver::initialize( coefficients );
}

void InverseSolver::solveImpl()
{
    LAMA_REGION( "Solver.Inverse.solve" );

    InverseSolverRuntime& runtime = getRuntime();

    LAMA_ASSERT_ERROR( runtime.mInverse.get(), "solve, but mInverse is NULL" );

    logStartSolve();
    *runtime.mSolution = ( *runtime.mInverse ) * ( *runtime.mRhs );
    logEndSolve();
}

void InverseSolver::setContext( ContextPtr context )
{
    Solver::setContext( context );
    getRuntime().mInverse->setContext( mContext );
}

void InverseSolver::computeInverse( Matrix& matrix ) const
{
    const IndexType n = matrix.getNumRows();
    boost::scoped_array<IndexType> permutation( new IndexType[n] );

    if ( typeid( matrix ) == typeid(DenseMatrix<float> ) )
    {
        DenseMatrix<float>& inverse = dynamic_cast<DenseMatrix<float>&>( matrix );

        invert( inverse, permutation.get() );
    }
    else if ( typeid( matrix ) == typeid(DenseMatrix<double> ) )
    {
        DenseMatrix<double>& inverse = dynamic_cast<DenseMatrix<double>&>( matrix );

        invert( inverse, permutation.get() );
    }
    else
    {
        //TODO: Implement fall back? (Create new Dense, invert new Dense and assign result to passed matrix?)
        LAMA_THROWEXCEPTION(
            "Computation of inverse is not supported for " << matrix << " because of a type missmatch (Type is not DenseMatrix)." );
    }
}

template<typename T>
void InverseSolver::decompose( DenseMatrix<T>& matrix, IndexType* const permutation ) const
{
    if ( matrix.getNumRows() != matrix.getNumColumns() )
    {
        LAMA_THROWEXCEPTION( "Can not decompose the not square matrix " << matrix );
    }
    if ( matrix.getDistribution().getNumPartitions() == 1 && matrix.getColDistribution().getNumPartitions() == 1 )
    {
        DenseStorage<T>& denseStorage = matrix.getLocalStorage();

        HostWriteAccess<T> denseValues( denseStorage.getData(), true );

        int error = lama_GETRF_cpu( CblasRowMajor, denseStorage.getNumRows(), denseStorage.getNumColumns(),
                                    denseValues.get(), denseStorage.getNumColumns(), &permutation[0] );
        if ( error != 0 )
        {
            LAMA_THROWEXCEPTION( "lama_GETRF_cpu failed" );
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Decomposition is not supported, because " << matrix << " is distributed." );
    }
}

template<typename T>
void InverseSolver::invert( DenseMatrix<T>& matrix, IndexType* const permutation ) const
{
    LAMA_REGION( "Solver.Inverse.invert" );

    typedef T ValueType;

    ContextPtr context = getCoefficients().getContextPtr();
    const LAMAInterface* lamaInterface = LAMAInterfaceRegistry::getRegistry().getInterface( context->getType() );

    if ( matrix.getNumRows() != matrix.getNumColumns() )
    {
        LAMA_THROWEXCEPTION( "Can not invert the not square matrix " << matrix );
    }
    if ( matrix.getDistribution().getNumPartitions() == 1 && matrix.getColDistribution().getNumPartitions() == 1 )
    {
        DenseStorage<T>& denseStorage = matrix.getLocalStorage();

        HostWriteAccess<T> denseValues( denseStorage.getData() );

        int error = lamaInterface->getLAPACKInterface<T>().getrf( CblasRowMajor, denseStorage.getNumRows(),
                    denseStorage.getNumColumns(), denseValues.get(),
                    denseStorage.getNumColumns(), &permutation[0] );

        if ( error != 0 )
        {
            LAMA_THROWEXCEPTION( "lama_GETRF_cpu failed" );
        }

        error = lamaInterface->getLAPACKInterface<T>().getri( CblasRowMajor, denseStorage.getNumRows(),
                denseValues.get(), denseStorage.getNumColumns(),
                &permutation[0] );

        if ( error != 0 )
        {
            LAMA_THROWEXCEPTION( "lama_GETRI_cpu failed" );
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Inversion is not supported, because " << matrix << " is distributed." );
    }
}

void InverseSolver::logStartSolve()
{
    mLogger->startTimer( "SolutionTimer" );
}

void InverseSolver::logEndSolve()
{
    L2Norm l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Final " );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );

}

InverseSolver::InverseSolverRuntime& InverseSolver::getRuntime()
{
    return mInverseSolverRuntime;
}

const InverseSolver::InverseSolverRuntime& InverseSolver::getConstRuntime() const
{
    return mInverseSolverRuntime;
}

SolverPtr InverseSolver::copy()
{
    return SolverPtr( new InverseSolver( *this ) );
}

}
