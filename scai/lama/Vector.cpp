/**
 * @file Vector.cpp
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
 * @brief Implementation of methods for the abstract class Vector.
 * @author Thomas Brandes
 * @date 31.10.2017
 */

#include <scai/lama/Vector.hpp>

#include <scai/common/TypeTraits.hpp>

namespace scai
{

using common::TypeTraits;

namespace lama
{

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::Vector( const IndexType size, hmemo::ContextPtr context ) :
       
   _Vector( size, context )
    
{
    // context will be set by base class _Vector

    SCAI_LOG_DEBUG( logger, "Vector<" << TypeTraits<ValueType>::id() 
                           << ">( size = " << size << ", ctx = " << getContext() << " )" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context ) :

    _Vector( distribution, context )

{
    // context will be set by base class _Vector

    SCAI_LOG_DEBUG( logger, "Vector<" << TypeTraits<ValueType>::id()
                           << ">( dist = " << getDistribution() << ", ctx = " << getContext() << " )" )
}

template<typename ValueType>
Vector<ValueType>::Vector( const _Vector& other ) : _Vector( other )
{
}

template<typename ValueType>
Vector<ValueType>::Vector( const Vector<ValueType>& other ) : _Vector( other )
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::~Vector()
{
    SCAI_LOG_DEBUG( logger, "~Vector<" << TypeTraits<ValueType>::id() << ">" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType Vector<ValueType>::getValueType() const
{
    return TypeTraits<ValueType>::stype;
}

} /* end namespace lama */

} /* end namespace scai */
