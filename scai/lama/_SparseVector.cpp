/**
 * @file _SparseVector.cpp
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
 * @brief Implementations of constructors/methods for class SparseVector.
 * @author Thomas Brandes
 * @date 16.01.2017
 */

// hpp
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/DenseVector.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/PartitionIO.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/common/BinaryOp.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unsupported.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>

// std
#include <ostream>

namespace scai
{

using common::TypeTraits;
using utilskernel::HArrayUtils;
using utilskernel::LArray;

using namespace hmemo;
using namespace dmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( _SparseVector::logger, "Vector.SparseVector" )

/* ------------------------------------------------------------------------- */
/*  Implementation of methods/constructors for _SparseVector                 */
/* ------------------------------------------------------------------------- */

_SparseVector::_SparseVector( const IndexType n ) :

    Vector( n )
{
}

_SparseVector::_SparseVector( const IndexType n, ContextPtr context ) :

    Vector( n, context )
{
}

_SparseVector::_SparseVector( DistributionPtr dist ) :

    Vector( dist )
{
}

_SparseVector::_SparseVector( DistributionPtr dist, ContextPtr context ) :

    Vector( dist, context )
{
}

_SparseVector::_SparseVector( const _SparseVector& other ) :

    Vector( other )
{
}

_SparseVector::_SparseVector( const Vector& other ) :

    Vector( other )
{
}

/* ------------------------------------------------------------------------- */
/*  Implementation of methods/constructors for _SparseVector                 */
/* ------------------------------------------------------------------------- */

_SparseVector* _SparseVector::create( common::scalar::ScalarType type )
{
    // There is only one factor for all vectors

    Vector* v = Vector::create( VectorCreateKeyType( Vector::SPARSE, type ) );

    // reinterpret cast is safe

    return reinterpret_cast<_SparseVector*>( v );
}

} /* end namespace lama */

} /* end namespace scai */
