/**
 * @file _DenseVector.cpp
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
 * @brief Implementations and instantiations for class DenseVector.
 * @author Thomas Brandes
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/_DenseVector.hpp>

// std
#include <ostream>

namespace scai
{

using namespace hmemo;
using namespace dmemo;

namespace lama
{

// common logger for all types of DenseVector

SCAI_LOG_DEF_LOGGER( _DenseVector::logger, "Vector.DenseVector" )

_DenseVector::_DenseVector( const IndexType n ) :

    Vector( n )
{
}

_DenseVector::_DenseVector( const IndexType n, ContextPtr context ) :

    Vector( n, context )
{
}

_DenseVector::_DenseVector( DistributionPtr dist ) :

    Vector( dist )
{
}

_DenseVector::_DenseVector( DistributionPtr dist, ContextPtr context ) :

    Vector( dist, context )
{
}

_DenseVector::_DenseVector( const _DenseVector& other ) :

    Vector( other )
{
}

_DenseVector::_DenseVector( const Vector& other ) :

    Vector( other )
{
}

_DenseVector* _DenseVector::create( common::ScalarType type )
{
    // There is only one factor for all vectors

    Vector* v = Vector::create( VectorCreateKeyType( Vector::DENSE, type ) );

    // reinterpret cast is safe

    return reinterpret_cast<_DenseVector*>( v );
}

} /* end namespace lama */

} /* end namespace scai */
