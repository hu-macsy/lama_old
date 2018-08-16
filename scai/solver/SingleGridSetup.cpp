/**
 * @file solver/SingleGridSetup.cpp
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
 * @brief Implementation of methods for the class SingleGridSetup.
 * @author Jiri Kraus
 * @date 27.10.2011
 */

// hpp
#include <scai/solver/SingleGridSetup.hpp>

// local library
#include <scai/solver/Jacobi.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

// tracing
#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SingleGridSetup<ValueType>::logger, "AMGSetup.SingleGridSetup" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_AMGSetup* SingleGridSetup<ValueType>::create()
{
    return (_AMGSetup*) ( new SingleGridSetup<ValueType>() );
}

template<typename ValueType>
AMGSetupCreateKeyType SingleGridSetup<ValueType>::createValue()
{
    return AMGSetupCreateKeyType( common::getScalarType<ValueType>(), "SingleGridSetup" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
SingleGridSetup<ValueType>::SingleGridSetup()
{
    SCAI_LOG_DEBUG( logger, "SingleGridSetup" )
}

template<typename ValueType>
SingleGridSetup<ValueType>::~SingleGridSetup()
{
}

/* ========================================================================= */
/*    Initialization                                                         */
/* ========================================================================= */

template<typename ValueType>
void SingleGridSetup<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "initialize_SingleGridSetup" )
    SCAI_LOG_DEBUG( logger, "SingleGridSetup::initialize" )

    // set default solver

    if ( !mSolver )
    {
        SCAI_LOG_DEBUG( logger, "new Jacobi" )
        Jacobi<ValueType>* jacobiSolver = new Jacobi<ValueType>( "10x SingleGridSetup Jacobi Solver" );
        CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 10 ) );
        jacobiSolver->setStoppingCriterion( criterion );
        mSolver.reset( jacobiSolver );
    }

    SCAI_LOG_DEBUG( logger, "mSolver->initialize" )
    mSolver->initialize( coefficients );
    SCAI_LOG_DEBUG( logger, "mIdentity.reset" )
    mIdentity.reset( coefficients.newMatrix() );
    SCAI_LOG_DEBUG( logger, "before identity" )
    mIdentity->setIdentity( coefficients.getRowDistributionPtr() );
    SCAI_LOG_DEBUG( logger, "after identity" )
    SCAI_LOG_DEBUG( logger, "Identity matrix = " << *mIdentity )

    mSolutionVector.reset( coefficients.newTargetVector() );
    mRhsVector.reset( coefficients.newTargetVector() );
    mTmpResVector.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    Getter methods                                                         */
/* ========================================================================= */

template<typename ValueType>
Solver<ValueType>& SingleGridSetup<ValueType>::getCoarseLevelSolver()
{
    return *mSolver;
}

template<typename ValueType>
IndexType SingleGridSetup<ValueType>::getNumLevels()
{
    return 2;
}

template<typename ValueType>
Solver<ValueType>& SingleGridSetup<ValueType>::getSmoother( const IndexType )
{
    return *mSolver;
}

template<typename ValueType>
const Matrix<ValueType>& SingleGridSetup<ValueType>::getGalerkin( const IndexType )
{
    return mSolver->getCoefficients();
}

template<typename ValueType>
const Matrix<ValueType>& SingleGridSetup<ValueType>::getRestriction( const IndexType )
{
    return *mIdentity;
}

template<typename ValueType>
const Matrix<ValueType>& SingleGridSetup<ValueType>::getInterpolation( const IndexType )
{
    return *mIdentity;
}

template<typename ValueType>
Vector<ValueType>& SingleGridSetup<ValueType>::getSolutionVector( const IndexType )
{
    return *mSolutionVector;
}

template<typename ValueType>
Vector<ValueType>& SingleGridSetup<ValueType>::getRhsVector( const IndexType )
{
    return *mRhsVector;
}

template<typename ValueType>
Vector<ValueType>& SingleGridSetup<ValueType>::getTmpResVector( const IndexType )
{
    return *mTmpResVector;
}

template<typename ValueType>
std::string SingleGridSetup<ValueType>::getCouplingPredicateInfo() const
{
    return "No coupling predicate.";
}

template<typename ValueType>
std::string SingleGridSetup<ValueType>::getColoringInfo() const
{
    return "No coloring.";
}

template<typename ValueType>
std::string SingleGridSetup<ValueType>::getInterpolationInfo() const
{
    return "Identity.";
}

template<typename ValueType>
std::string SingleGridSetup<ValueType>::getSmootherInfo() const
{
    return mSolver->getId();
}

template<typename ValueType>
std::string SingleGridSetup<ValueType>::getCoarseLevelSolverInfo() const
{
    return mSolver->getId();
}

template<typename ValueType>
void SingleGridSetup<ValueType>::setCoarseLevelSolver( SolverPtr<ValueType> solver )
{
    mSolver = solver;
}

template<typename ValueType>
void SingleGridSetup<ValueType>::setSmoother( SolverPtr<ValueType> solver )
{
    mSolver = solver;
}

template<typename ValueType>
void SingleGridSetup<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "SingleGridSetup";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SingleGridSetup, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
