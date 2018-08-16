/**
 * @file common/safer_memcpy.hpp
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
 * @brief A safe alternative to memcpy when size is zero
 * @author Andreas Longva
 * @date 20.10.2017
 */

#pragma once

#include <algorithm>

#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace common
{
    /** A safe alternative to memcpy when the size of the data is zero.
     * 
     * According to the C++ standard, std::memcpy(dest, src, count) is
     * undefined behavior if dest or src is a null pointer, even when
     * count is zero.
     * 
     * safer_memcpy is a safe alternative to memcpy, for which
     * it is safe to pass null pointers when count is zero
     * (but otherwise not).
     */
    inline void * safer_memcpy( void * dest, const void * src, std::size_t count )
    {
        SCAI_ASSERT_DEBUG( count == 0 || (dest != NULL && src != NULL),
                           "Null pointers are only allowed when count is zero.");
        
        // We cannot use std::copy with void * directly, but by casting the void pointer
        // to a char pointer, we can.
        const char * begin = static_cast<const char *>( src );
        std::copy(begin, begin + count, static_cast<char *>( dest ) );
        return dest;
    }
} // common

} // scai
