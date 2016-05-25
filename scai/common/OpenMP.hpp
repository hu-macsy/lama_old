/**
 * @file OpenMP.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Common defintions for optional use of OpenMP
 * @author Thomas Brandes
 * @date 11.06.2013
 */

#pragma once

#ifdef _OPENMP

#include <omp.h>

#else

#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_set_num_threads( x ) 

#endif

/** atomicAdd used for reductions as reduction directive is unsupported for complex numbers.
 *
 *  Note: template specialization used for float and double
 */

template<typename ValueType>
inline void atomicAdd( ValueType& sharedResult, const ValueType& threadResult )
{
    #pragma omp critical
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( float& sharedResult, const float& threadResult )
{
    #pragma omp atomic
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( double& sharedResult, const double& threadResult )
{
    #pragma omp atomic
    sharedResult += threadResult;
}
