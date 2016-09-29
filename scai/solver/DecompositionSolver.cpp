/**
 * @file DecompositionSolver.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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

// std
#include <sstream>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( DecompositionSolver::logger, "Solver.DecompositionSolver" )

DecompositionSolver::DecompositionSolver( const std::string& id )
    : Solver( id )
{
    SCAI_LOG_INFO( DecompositionSolver::logger, "DecompositionSolver, id = " << id )
}

DecompositionSolver::DecompositionSolver( const std::string& id, LoggerPtr logger )
    : Solver( id, logger )
{
    SCAI_LOG_INFO( DecompositionSolver::logger, "DecompositionSolver, id = " << id )
}

DecompositionSolver::DecompositionSolver( const DecompositionSolver& other )
    : Solver( other )
{
    SCAI_LOG_INFO( DecompositionSolver::logger, "DecompositionSolver, id = " << other.mId )
}

DecompositionSolver::DecompositionSolverRuntime::DecompositionSolverRuntime()
    : SolverRuntime()
{
}

DecompositionSolver::~DecompositionSolver()
{
    SCAI_LOG_INFO( logger, "~DecompositionSolver" )
}

DecompositionSolver::DecompositionSolverRuntime::~DecompositionSolverRuntime()
{
}

/* --------------------------------------------------------------------------- */

void DecompositionSolver::initialize( const lama::Matrix& coefficients )
{
    SCAI_REGION( "Solver.DecompositionSolver.intialize" )
    SCAI_LOG_INFO( logger, "Initializing with " << coefficients )
    Solver::initialize( coefficients );
    // TODO: check symmetry
    DecompositionSolverRuntime& runtime = getRuntime();
    runtime.mIsSymmetric = false;
}

/* --------------------------------------------------------------------------- */

void DecompositionSolver::solveImpl()
{
	SCAI_REGION( "Solver.DecompositionSolver.solveImpl" )
    // for each supported arithmetic type we have to dynamic cast and instantiate typed version
#define SCAI_SOLVER_TYPE_CAST( _type )                                                                                  \
    {                                                                                               \
        const lama::SparseMatrix<_type>* sparseTypedCoefficients =                                  \
                dynamic_cast<const lama::SparseMatrix<_type>*>( getRuntime().mCoefficients );       \
        if ( sparseTypedCoefficients )                                                              \
        {                                                                                           \
            solveImplTyped( *sparseTypedCoefficients );                                             \
            return;                                                                                 \
        }                                                                                           \
    }
    SCAI_COMMON_LOOP( SCAI_SOLVER_TYPE_CAST, SCAI_NUMERIC_TYPES_HOST )
#undef SCAI_SOLVER_TYPE_CAST
    // has already been check in initialize, but in any case
    COMMON_THROWEXCEPTION        (
        getConstRuntime().mCoefficients << ": unsupported matrix type (only SparseMatrix<ValueType> supported)." )
}

/* --------------------------------------------------------------------------- */

template <typename ValueType>
void DecompositionSolver::solveImplTyped( const lama::SparseMatrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.DecompositionSolver.solve" )
    SCAI_LOG_INFO( logger, "solveImplTyped for ValueType=" << scai::common::TypeTraits<ValueType>::stype )

    DecompositionSolverRuntime& runtime = getRuntime();
    logStartSolve();

    // only implemented for replicated matrices

    if( coefficients.getRowDistribution().isReplicated() && coefficients.getColDistribution().isReplicated() )
    {
    	// need CSR storage
    	const lama::CSRStorage<ValueType>* csrStorage;

    	if( coefficients.getLocalStorage().getFormat() == lama::Format::CSR )
    	{
    		csrStorage = dynamic_cast<const lama::CSRStorage<ValueType>* >( &coefficients.getLocalStorage() );
    	}
    	else
    	{
    		csrStorage = new lama::CSRStorage<ValueType>( coefficients.getLocalStorage() );
    	}

        lama::DenseVector<ValueType>& denseSolution = dynamic_cast<lama::DenseVector<ValueType>&>( *getRuntime().mSolution );
        const lama::DenseVector<ValueType>& denseRhs = dynamic_cast<const lama::DenseVector<ValueType>&>( *getRuntime().mRhs );

    	// call decompostion
    	IndexType numRows = csrStorage->getNumRows();
    	IndexType nnz = csrStorage->getNumValues();

	    hmemo::ContextPtr prefContext = csrStorage->getContextPtr();
    	kregistry::KernelTraitContextFunction<sparsekernel::CSRKernelTrait::decomposition<ValueType> > decomposition;
    	hmemo::ContextPtr loc = hmemo::Context::getContextPtr( decomposition.validContext( prefContext->getType() ) );

    	std::cout << "prefContext=" << *prefContext << " loc=" << *loc << std::endl;

    	hmemo::ReadAccess<IndexType> rCSRIA( csrStorage->getIA(), loc );
        hmemo::ReadAccess<IndexType> rCSRJA( csrStorage->getJA(), loc );
        hmemo::ReadAccess<ValueType> rCSRValues( csrStorage->getValues(), loc );
        hmemo::ReadAccess<ValueType> rRHS( denseRhs.getLocalValues(), loc );
        hmemo::WriteOnlyAccess<ValueType> wSol( denseSolution.getLocalValues(), loc, numRows );

        SCAI_CONTEXT_ACCESS( loc );
        decomposition[loc->getType()]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
            rRHS.get(), numRows, nnz, runtime.mIsSymmetric );
    }
    else
    {
    	COMMON_THROWEXCEPTION( "DecompositionSolver not implemented for distributed matrices yet." )
    }

    logEndSolve();
}

/* --------------------------------------------------------------------------- */

void DecompositionSolver::logStartSolve()
{
    mLogger->startTimer( "SolutionTimer" );
}

/* --------------------------------------------------------------------------- */

void DecompositionSolver::logEndSolve()
{
    lama::L2Norm l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Final " );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );
}

/* --------------------------------------------------------------------------- */

DecompositionSolver::DecompositionSolverRuntime& DecompositionSolver::getRuntime()
{
    return mDecompositionSolverRuntime;
}

/* --------------------------------------------------------------------------- */

const DecompositionSolver::DecompositionSolverRuntime& DecompositionSolver::getConstRuntime() const
{
    return mDecompositionSolverRuntime;
}

/* --------------------------------------------------------------------------- */

SolverPtr DecompositionSolver::copy()
{
    return SolverPtr( new DecompositionSolver( *this ) );
}

/* --------------------------------------------------------------------------- */

void DecompositionSolver::writeAt( std::ostream& stream ) const
{
    stream << "DecompositionSolver ( id = " << mId << " )";
}

/* --------------------------------------------------------------------------- */

std::string DecompositionSolver::createValue()
{
    return "DecompositionSolver";
}

Solver* DecompositionSolver::create( const std::string name )
{
    return new DecompositionSolver( name );
}

} /* end namespace solver */

} /* end namespace scai */
