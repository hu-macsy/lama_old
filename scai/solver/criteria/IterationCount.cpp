/**
 * @file IterationCount.cpp
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
 * @brief IterationCount.cpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

// hpp
#include <scai/solver/criteria/IterationCount.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, IterationCount<ValueType>::logger, "Criterion.IterationCount" );

template<typename ValueType>
IterationCount<ValueType>::IterationCount() : 

    Criterion<ValueType>(), 
    mIterationExtrema( 1 )
{
}

template<typename ValueType>
IterationCount<ValueType>::IterationCount( const IndexType iterationExtrema ) : 

    Criterion<ValueType>(), 
    mIterationExtrema( iterationExtrema )
{
    SCAI_LOG_DEBUG( logger, "Creating IterationCount with " << iterationExtrema );
}

template<typename ValueType>
IterationCount<ValueType>::IterationCount( const IterationCount<ValueType>& other ) : 

    Criterion<ValueType>( other ), 
    mIterationExtrema( other.mIterationExtrema )
{
}

template<typename ValueType>
IterationCount<ValueType>::~IterationCount()
{
}

template<typename ValueType>
bool IterationCount<ValueType>::isSatisfied( const IterativeSolver<ValueType>& solver )
{
    SCAI_LOG_INFO( logger,
                   "Iteration Extrema = " << mIterationExtrema << ", Iteration Count = " << solver.getIterationCount() );
    return solver.getIterationCount() >= mIterationExtrema;
}

template<typename ValueType>
IndexType IterationCount<ValueType>::getIterationExtrema() const
{
    return mIterationExtrema;
}

template<typename ValueType>
void IterationCount<ValueType>::setIterationExtrema( IndexType iterationExtrema )
{
    mIterationExtrema = iterationExtrema;
}

template<typename ValueType>
void IterationCount<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ItCount<" << getIterationExtrema() << ">";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( IterationCount, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
