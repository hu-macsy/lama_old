/**
 * @file SolutionProxy.cpp
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
 * @brief SolutionProxy.cpp
 * @author Jiri Kraus
 * @date 07.06.2011
 */

// hpp
#include <scai/solver/SolutionProxy.hpp>

namespace scai
{

namespace solver
{

SolutionProxy::SolutionProxy()
    : mIsDirty( true )
{
}

SolutionProxy::SolutionProxy( lama::Vector* const solution )
    : mSolution( solution ), mIsDirty( true )
{
}

SolutionProxy::~SolutionProxy()
{
}

const lama::Vector& SolutionProxy::getConstReference() const
{
    return ( *mSolution );
}

lama::Vector& SolutionProxy::operator*()
{
    return getReference();
}

void SolutionProxy::operator=( lama::Vector* const newVector )
{
    setDirty( true );
    mSolution = newVector;
}

bool SolutionProxy::isDirty() const
{
    return mIsDirty;
}

void SolutionProxy::setDirty( bool isDirty )
{
    mIsDirty = isDirty;
}

lama::Vector& SolutionProxy::getReference()
{
    setDirty( true );
    return *mSolution;
}

} /* end namespace solver */

} /* end namespace scai */
