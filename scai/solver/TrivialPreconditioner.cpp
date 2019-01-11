/**
 * @file TrivialPreconditioner.cpp
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
 * @brief Implementation of methods for the class TrivialPreconditioner.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/TrivialPreconditioner.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using lama::Matrix;
using lama::Vector;

namespace solver
{

template<typename ValueType>
TrivialPreconditioner<ValueType>::TrivialPreconditioner( const std::string& id )
    : Solver<ValueType>( id )
{
}

template<typename ValueType>
TrivialPreconditioner<ValueType>::TrivialPreconditioner( const std::string& id, LoggerPtr logger )
    : Solver<ValueType>( id, logger )
{
}

template<typename ValueType>
TrivialPreconditioner<ValueType>::TrivialPreconditioner( const TrivialPreconditioner& other )
    : Solver<ValueType>( other )
{
}

template<typename ValueType>
TrivialPreconditioner<ValueType>::TrivialPreconditionerRuntime::TrivialPreconditionerRuntime()
{
}

template<typename ValueType>
TrivialPreconditioner<ValueType>::~TrivialPreconditioner()
{
}

template<typename ValueType>
TrivialPreconditioner<ValueType>::TrivialPreconditionerRuntime::~TrivialPreconditionerRuntime()
{
}

template<typename ValueType>
void TrivialPreconditioner<ValueType>::solveImpl()
{
    lama::Vector<ValueType>& solution  = getRuntime().mSolution.getReference();
    const lama::Vector<ValueType>& rhs = *getRuntime().mRhs;

    solution = rhs;
}

template<typename ValueType>
TrivialPreconditioner<ValueType>* TrivialPreconditioner<ValueType>::copy()
{
    return new TrivialPreconditioner( *this );
}

template<typename ValueType>
typename TrivialPreconditioner<ValueType>::TrivialPreconditionerRuntime& TrivialPreconditioner<ValueType>::getRuntime()
{
    return mTrivialPreconditionerRuntime;
}

template<typename ValueType>
const typename TrivialPreconditioner<ValueType>::TrivialPreconditionerRuntime& TrivialPreconditioner<ValueType>::getRuntime() const
{
    return mTrivialPreconditionerRuntime;
}

template<typename ValueType>
std::string TrivialPreconditioner<ValueType>::createValue()
{
    return "TrivialPreconditioner";
}

template<typename ValueType>
void TrivialPreconditioner<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "TrivialPreconditioner ( id = " << this->getId() << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( TrivialPreconditioner, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
