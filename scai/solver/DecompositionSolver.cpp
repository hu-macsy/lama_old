/**
 * @file DecompositionSolver.cpp
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
 * @brief DecompositionSolver.cpp
 * @author Lauretta Schubert
 * @date 20.07.2016
 */

// hpp
#include <scai/solver/DecompositionSolver.hpp>

// internal scai libraries
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/kregistry.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <sstream>
#include <memory>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, 
                              DecompositionSolver<ValueType>::logger, "Solver.DecompositionSolver" )

using lama::Matrix;
using lama::Vector;
using lama::DenseVector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* DecompositionSolver<ValueType>::create()
{
    return new DecompositionSolver<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType DecompositionSolver<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "DecompositionSolver" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
DecompositionSolver<ValueType>::DecompositionSolver( const std::string& id ) : 

    Solver<ValueType>( id )

{
    SCAI_LOG_INFO( DecompositionSolver<ValueType>::logger, "DecompositionSolver, id = " << id )
}

template<typename ValueType>
DecompositionSolver<ValueType>::DecompositionSolver( const std::string& id, LoggerPtr logger ) : 

    Solver<ValueType>( id, logger )
{
    SCAI_LOG_INFO( DecompositionSolver<ValueType>::logger, "DecompositionSolver, id = " << id )
}

template<typename ValueType>
DecompositionSolver<ValueType>::DecompositionSolver( const DecompositionSolver<ValueType>& other ) : 
 
    Solver<ValueType>( other )
{
    SCAI_LOG_INFO( DecompositionSolver<ValueType>::logger, "DecompositionSolver, id = " << other.getId() )
}

template<typename ValueType>
DecompositionSolver<ValueType>::~DecompositionSolver()
{
    SCAI_LOG_INFO( logger, "~DecompositionSolver" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void DecompositionSolver<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.DecompositionSolver.intialize" )

    SCAI_LOG_INFO( logger, "Initializing with " << coefficients )

    Solver<ValueType>::initialize( coefficients );

    // TODO: check symmetry

    DecompositionSolverRuntime& runtime = getRuntime();
    runtime.mIsSymmetric = false;
}

/* ========================================================================= */
/*    solve : implemantation                                                 */
/* ========================================================================= */

template <typename ValueType>
void DecompositionSolver<ValueType>::solveImpl()
{
    SCAI_REGION( "Solver.DecompositionSolver.solve" )

    DecompositionSolverRuntime& runtime = getRuntime();

    logStartSolve();

    const Matrix<ValueType>& coefficients = *runtime.mCoefficients;

    // only implemented for replicated matrices

    if ( !coefficients.getRowDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "DecompositionSolver not implemented for distributed matrices yet." )
    }

    if ( !coefficients.getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "DecompositionSolver not implemented for matrices with distributed columns yet." )
    }

    // from now on we do it all on storage, but we need CSR storage for the solver

    std::unique_ptr<lama::CSRStorage<ValueType> > tmpCSRStorage;

    const lama::CSRStorage<ValueType>* csrStorage;

    if ( coefficients.getLocalStorage().getFormat() == lama::Format::CSR )
    {
        csrStorage = dynamic_cast<const lama::CSRStorage<ValueType>* >( &coefficients.getLocalStorage() );
    }
    else
    {
        tmpCSRStorage.reset( new lama::CSRStorage<ValueType>() );
        tmpCSRStorage->assign( coefficients.getLocalStorage() );
        SCAI_LOG_WARN( logger, "new tmp csr storage = " << *tmpCSRStorage )
        csrStorage = tmpCSRStorage.get();
    }

    // we must sort the column indexes for solving the matrix, no diagonal flag

    lama::CSRStorage<ValueType>* xCSRStorage = const_cast<lama::CSRStorage<ValueType>*>( csrStorage );
    xCSRStorage->sortRows();

    SCAI_LOG_INFO( logger, "csrStorage = " << *csrStorage )

    Vector<ValueType>& solution = getRuntime().mSolution.getReference();  // -> dirty
    const Vector<ValueType>& rhs = *getRuntime().mRhs;

    DenseVector<ValueType>& denseSolution = reinterpret_cast<DenseVector<ValueType>&>( solution );
    const DenseVector<ValueType>& denseRhs = reinterpret_cast<const lama::DenseVector<ValueType>&>( rhs );

    SCAI_LOG_INFO( logger, "solution = " << denseSolution << ", rhs = " << denseRhs )

    // call decompostion

    IndexType numRows = csrStorage->getNumRows();
    IndexType nnz = csrStorage->getNumValues();

    hmemo::ContextPtr prefContext = csrStorage->getContextPtr();

    kregistry::KernelTraitContextFunction<sparsekernel::CSRKernelTrait::decomposition<ValueType> > decomposition;

    hmemo::ContextPtr loc = hmemo::Context::getContextPtr( decomposition.validContext( prefContext->getType() ) );

    SCAI_LOG_INFO( logger, "prefContext=" << *prefContext << " loc=" << *loc )

    {
        hmemo::ReadAccess<IndexType> rCSRIA( csrStorage->getIA(), loc );
        hmemo::ReadAccess<IndexType> rCSRJA( csrStorage->getJA(), loc );
        hmemo::ReadAccess<ValueType> rCSRValues( csrStorage->getValues(), loc );
        hmemo::ReadAccess<ValueType> rRHS( denseRhs.getLocalValues(), loc );

        hmemo::WriteOnlyAccess<ValueType> wSol( denseSolution.getLocalValues(), loc, numRows );

        SCAI_CONTEXT_ACCESS( loc );

        decomposition[loc->getType()]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
                                       rRHS.get(), numRows, nnz, runtime.mIsSymmetric );
    }

    // Note: accesses must be released before the next call

    logEndSolve();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DecompositionSolver<ValueType>::logStartSolve()
{
    mLogger->startTimer( "SolutionTimer" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DecompositionSolver<ValueType>::logEndSolve()
{
    lama::L2Norm<ValueType> l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Final " );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );
}

/* ========================================================================= */
/*       Runtime                                                             */
/* ========================================================================= */

template<typename ValueType>
typename DecompositionSolver<ValueType>::DecompositionSolverRuntime& DecompositionSolver<ValueType>::getRuntime()
{
    return mDecompositionSolverRuntime;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const typename DecompositionSolver<ValueType>::DecompositionSolverRuntime& DecompositionSolver<ValueType>::getRuntime() const
{
    return mDecompositionSolverRuntime;
}

/* ========================================================================= */
/*       Virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
DecompositionSolver<ValueType>* DecompositionSolver<ValueType>::copy()
{
    return new DecompositionSolver( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DecompositionSolver<ValueType>::writeAt( std::ostream& stream ) const
{
    const char* typeId = common::TypeTraits<ValueType>::id();

    stream << "DecompositionSolver<" << typeId << "> ( id = " << this->getId() << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DecompositionSolver, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
