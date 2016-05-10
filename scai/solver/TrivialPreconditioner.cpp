/**
 * @file TrivialPreconditioner.cpp
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
 * @brief TrivialPreconditioner.cpp
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/TrivialPreconditioner.hpp>

namespace scai
{

using lama::Matrix;

namespace solver
{

TrivialPreconditioner::TrivialPreconditioner( const std::string& id )
    : Solver( id )
{
}

TrivialPreconditioner::TrivialPreconditioner( const std::string& id, LoggerPtr logger )
    : Solver( id, logger )
{
}

TrivialPreconditioner::TrivialPreconditioner( const TrivialPreconditioner& other )
    : Solver( other )
{
}

TrivialPreconditioner::TrivialPreconditionerRuntime::TrivialPreconditionerRuntime()
    : SolverRuntime()
{
}

TrivialPreconditioner::~TrivialPreconditioner()
{
}

TrivialPreconditioner::TrivialPreconditionerRuntime::~TrivialPreconditionerRuntime()
{
}

void TrivialPreconditioner::initialize( const Matrix& coefficients )
{
    Solver::initialize( coefficients );
}

void TrivialPreconditioner::solveImpl()
{
    *( getRuntime().mSolution ) = *( getRuntime().mRhs );
}

SolverPtr TrivialPreconditioner::copy()
{
    return SolverPtr( new TrivialPreconditioner( *this ) );
}

TrivialPreconditioner::TrivialPreconditionerRuntime& TrivialPreconditioner::getRuntime()
{
    return mTrivialPreconditionerRuntime;
}

const TrivialPreconditioner::TrivialPreconditionerRuntime& TrivialPreconditioner::getConstRuntime() const
{
    return mTrivialPreconditionerRuntime;
}

std::string TrivialPreconditioner::createValue()
{
	return "TrivialPreconditioner";
}

Solver* TrivialPreconditioner::create( const std::string name )
{
	return new TrivialPreconditioner( name );
}

void TrivialPreconditioner::writeAt( std::ostream& stream ) const
{
    stream << "TrivialPreconditioner ( id = " << mId << " )";
}

} /* end namespace solver */

} /* end namespace scai */
