/**
 * @file SolutionProxy.cpp
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
 * @brief Implementation of methods for the class SolutionProxy.
 * @author Jiri Kraus
 * @date 07.06.2011
 */

#include <scai/solver/SolutionProxy.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType>
SolutionProxy<ValueType>::SolutionProxy() : mIsDirty( true )
{
}

template<typename ValueType>
SolutionProxy<ValueType>::SolutionProxy( lama::Vector<ValueType>* solution ) : 

    mSolution( solution ), 
    mIsDirty( true )
{
}

template<typename ValueType>
SolutionProxy<ValueType>::~SolutionProxy()
{
}

template<typename ValueType>
const lama::Vector<ValueType>& SolutionProxy<ValueType>::getConstReference() const
{
    return *mSolution;
}

template<typename ValueType>
void SolutionProxy<ValueType>::operator=( lama::Vector<ValueType>* newVector )
{
    setDirty( true );
    mSolution = newVector;
}

template<typename ValueType>
bool SolutionProxy<ValueType>::isDirty() const
{
    return mIsDirty;
}

template<typename ValueType>
void SolutionProxy<ValueType>::setDirty( bool isDirty )
{
    mIsDirty = isDirty;
}

template<typename ValueType>
lama::Vector<ValueType>& SolutionProxy<ValueType>::getReference()
{
    this->setDirty( true );
    return *mSolution;
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SolutionProxy, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
