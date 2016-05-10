/**
 * @file IterationCount.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief IterationCount.cpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

// hpp
#include <scai/solver/criteria/IterationCount.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( IterationCount::logger, "Criterion.IterationCount" );

IterationCount::IterationCount()
    : Criterion(), mIterationExtrema( 1 )
{
}

IterationCount::IterationCount( const IndexType iterationExtrema )
    : Criterion(), mIterationExtrema( iterationExtrema )
{
    SCAI_LOG_DEBUG( logger, "Creating IterationCount with " << iterationExtrema );
}

IterationCount::IterationCount( const IterationCount &other )
    : Criterion( other ), mIterationExtrema( other.mIterationExtrema )

{
}

IterationCount::~IterationCount()
{
}

bool IterationCount::isSatisfied( const IterativeSolver& solver )
{
    SCAI_LOG_INFO( logger,
                   "Iteration Extrema = " << mIterationExtrema << ", Iteration Count = " << solver.getIterationCount() );

    return solver.getIterationCount() >= mIterationExtrema;
}

IndexType IterationCount::getIterationExtrema() const
{
    return mIterationExtrema;
}

void IterationCount::setIterationExtrema( IndexType iterationExtrema )
{
    mIterationExtrema = iterationExtrema;
}

void IterationCount::writeAt( std::ostream& stream ) const
{
    stream << "ItCount<" << getIterationExtrema() << ">";
}

} /* end namespace solver */

} /* end namespace scai */
