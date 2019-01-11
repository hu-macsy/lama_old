/**
 * @file utilskernel/freeFunction.hpp
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
 * @brief Definiton of free functions to construct heterogeneous arrays
 * @author Thomas Brandes
 * @date 09.03.18
 */
#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

namespace scai
{

namespace utilskernel
{

template<typename ValueType> 
hmemo::HArray<ValueType> randomHArray( const IndexType size, const IndexType bound, hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    hmemo::HArray<ValueType> array( size );
    HArrayUtils::fillRandom( array, bound, 1.0f, ctx );
    return array;
}

template<typename ValueType> 
hmemo::HArray<ValueType> sparseRandomHArray( 
    const IndexType size, 
    const ValueType zero,
    const float fillRate,
    const IndexType bound, 
    hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    hmemo::HArray<ValueType> array;
    HArrayUtils::setSameValue( array, size, zero, ctx );
    HArrayUtils::fillRandom( array, bound, fillRate, ctx );
    return array;
}

template<typename ValueType> 
hmemo::HArray<ValueType> fillHArray( const IndexType size, const ValueType val, hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    hmemo::HArray<ValueType> array;
    HArrayUtils::setSameValue( array, size, val, ctx );
    return array;
}

template<typename ValueType> 
hmemo::HArray<ValueType> convertHArray( const hmemo::_HArray& other, hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    hmemo::HArray<ValueType> array;
    HArrayUtils::_assign( array, other, ctx );
    return array;
}

} /* end namespace utilskernel */

} /* end namespace scai */

