/**
 * @file common/mepr/TypeList.hpp
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
 * @brief TypeList structure
 * @author Eric Schricker
 * @date 03.03.2016
 */

#pragma once

#include <scai/common/macros/count.hpp>

namespace scai
{

namespace common
{

/** Own namespace for all stuff belonging to metaprogramming */

namespace mepr
{

/*
 * NullType, used for termination
 */
class NullType {};

/**
 * Definition of struct to build list of types.
 *
 * @tparam H head type
 * @tparam T struct of tail of the list, terminates with NullType
 */

template<typename H, typename T>
struct TypeList
{
    typedef H head;
    typedef T tail;
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */

/*
 * Internal used macros for TypeList creation
 *
 * TYPELIST( N, ... ) for usage. N is the number of types
 */

#define TYPELIST_0( ) scai::common::mepr::NullType
#define TYPELIST_1( T1 ) scai::common::mepr::TypeList<T1,scai::common::mepr::NullType>
#define TYPELIST_2( T1, T2 ) scai::common::mepr::TypeList<T1,TYPELIST_1( T2 ) >
#define TYPELIST_3( T1, T2, T3 ) scai::common::mepr::TypeList<T1,TYPELIST_2( T2, T3 ) >
#define TYPELIST_4( T1, T2, T3, T4 ) scai::common::mepr::TypeList<T1,TYPELIST_3( T2, T3, T4 ) >
#define TYPELIST_5( T1, T2, T3, T4, T5 ) scai::common::mepr::TypeList<T1,TYPELIST_4( T2, T3, T4, T5 ) >
#define TYPELIST_6( T1, T2, T3, T4, T5, T6 ) scai::common::mepr::TypeList<T1,TYPELIST_5( T2, T3, T4, T5, T6 ) >
#define TYPELIST_7( T1, T2, T3, T4, T5, T6, T7 ) scai::common::mepr::TypeList<T1,TYPELIST_6( T2, T3, T4, T5, T6, T7 ) >
#define TYPELIST_8( T1, T2, T3, T4, T5, T6, T7, T8 ) scai::common::mepr::TypeList<T1,TYPELIST_7( T2, T3, T4, T5, T6, T7, T8 ) >
#define TYPELIST_9( T1, T2, T3, T4, T5, T6, T7, T8, T9 ) scai::common::mepr::TypeList<T1,TYPELIST_8( T2, T3, T4, T5, T6, T7, T8, T9 ) >
#define TYPELIST_10( T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 ) scai::common::mepr::TypeList<T1,TYPELIST_9( T2, T3, T4, T5, T6, T7, T8, T9, T10 ) >
#define TYPELIST_11( T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11 ) scai::common::mepr::TypeList<T1,TYPELIST_10( T2, T3, T4, T5, T6, T7, T8, T9, T10, T11 ) >

// Help macro is required to guarantee that argument N is no more a macro itself

#define TYPELIST_N( N, ... ) TYPELIST_##N( __VA_ARGS__ )

/** Macro for type lists where the first argument N specifies the number of types.
 *  Be careful: the total number of arguments must be N + 1
 */

#define TYPELIST( N, ... ) TYPELIST_N( N, __VA_ARGS__ )

#define SCAI_TYPELIST( ... ) TYPELIST( SCAI_COMMON_COUNT_NARG( __VA_ARGS__ ), __VA_ARGS__ )
