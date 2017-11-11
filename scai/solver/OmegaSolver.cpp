/**
 * @file solver/OmegaSolver.cpp
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
 * @brief Test methods for the generic solver class OmegaSolver.
 * @author Kai Buschulte
 * @date 10.08.2011
 */

// hpp
#include <scai/solver/OmegaSolver.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( OmegaSolver::logger, "Solver.IterativeSolver.OmegaSolver" )

OmegaSolver::OmegaSolver( const std::string& id )
    : IterativeSolver( id ), mOmega( 0.5 )
{
}

OmegaSolver::OmegaSolver( const std::string& id, const lama::Scalar omega )
    : IterativeSolver( id ), mOmega( omega )
{
}

OmegaSolver::OmegaSolver( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger ), mOmega( 0.5 )
{
}

OmegaSolver::OmegaSolver( const std::string& id, const lama::Scalar omega, LoggerPtr logger )
    : IterativeSolver( id, logger ), mOmega( omega )
{
}

OmegaSolver::OmegaSolver( const OmegaSolver& other )
    : IterativeSolver( other ), mOmega( other.mOmega )
{
}

OmegaSolver::OmegaSolverRuntime::OmegaSolverRuntime()
    : IterativeSolverRuntime()
{
}

OmegaSolver::~OmegaSolver()
{
}

OmegaSolver::OmegaSolverRuntime::~OmegaSolverRuntime()
{
}

void OmegaSolver::setOmega( const ValueType omega )
{
    mOmega = omega;
}

ValueType OmegaSolver::getOmega() const
{
    return mOmega;
}

} /* end namespace solver */

} /* end namespace scai */
