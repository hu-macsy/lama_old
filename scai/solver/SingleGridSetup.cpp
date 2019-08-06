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
    SCAI_LOG_DEBUG( logger, "~SingleGridSetup" )
}

/* ========================================================================= */
/*    createSolver()                                                         */
/* ========================================================================= */

template<typename ValueType>
SolverPtr<ValueType> SingleGridSetup<ValueType>::createSolver( bool )
{
    // Note: we do not have a coarse level at all here, so isCoarseLevel does not matter

    auto jacobiSolver = std::make_shared<Jacobi<ValueType>>( "10x SingleGridSetup Jacobi Solver" );

    auto criterion = std::make_shared<IterationCount<ValueType>>( 10 );

    jacobiSolver->setStoppingCriterion( criterion );

    return jacobiSolver;
}

/* ========================================================================= */
/*    createMatrixHierarchy()                                                */
/* ========================================================================= */

template<typename ValueType>
void SingleGridSetup<ValueType>::createMatrixHierarchy()
{
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
void SingleGridSetup<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "SingleGridSetup( #levels = " << this->getNumLevels() << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SingleGridSetup, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
