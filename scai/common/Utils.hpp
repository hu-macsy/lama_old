/**
 * @file Utils.hpp
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
 * @brief Definition of struct that contains util routines overloaded for all possible index types
 * @author Thomas Brandes
 * @date 15.09.2016
 */

#pragma once

#include <scai/common/cuda/CUDACallable.hpp>

namespace scai
{

namespace common
{

/** Structure (instead of namespace) that contains the required utility functions
 *  for each possible index type.
 */

struct Utils
{
    /*
     * valid index verifies 0 <= index < size for signed types and index < size for unsigned types
     *
     * Reasoning: common code causes compiler warnings as comparison >= 0 is useless for unsigned types
     */
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const short& index, const short& size );
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const unsigned short& index, const unsigned short& size );
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const int& index, const int& size );
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const unsigned int& index, const unsigned int& size );
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const long& index, const long& size );
    static inline CUDA_CALLABLE_MEMBER bool validIndex( const unsigned long& index, const unsigned long& size );
};

// -------------------------------- validIndex -----------------------

bool Utils::validIndex( const short& index, const short& size )
{
    return index >= 0 && index < size;
}

bool Utils::validIndex( const unsigned short& index, const unsigned short& size )
{
    return index < size;
}

bool Utils::validIndex( const int& index, const int& size )
{
    return index >= 0 && index < size;
}

bool Utils::validIndex( const unsigned int& index, const unsigned int& size )
{
    return index < size;
}

bool Utils::validIndex( const long& index, const long& size )
{
    return index >= 0 && index < size;
}

bool Utils::validIndex( const unsigned long& index, const unsigned long& size )
{
    return index < size;
}

} /* end namespace common */

} /* end namespace scai */
