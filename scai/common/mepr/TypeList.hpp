/**
 * @file common/mepr/TypeList.hpp
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
 * @brief TypeList structure
 * @author Eric Schricker
 * @date 03.03.2016
 */

#pragma once

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
 * TYPELIST( NR, ... ) for usage. NR is the number of types
 */

#define TYPELIST_1( T1 ) scai::common::mepr::TypeList<T1,common::mepr::NullType>
#define TYPELIST_2( T1, T2 ) scai::common::mepr::TypeList<T1,TYPELIST_1( T2 ) >
#define TYPELIST_3( T1, T2, T3 ) scai::common::mepr::TypeList<T1,TYPELIST_2( T2, T3 ) >
#define TYPELIST_4( T1, T2, T3, T4 ) scai::common::mepr::TypeList<T1,TYPELIST_3( T2, T3, T4 ) >
#define TYPELIST_5( T1, T2, T3, T4, T5 ) scai::common::mepr::TypeList<T1,TYPELIST_4( T2, T3, T4, T5 ) >
#define TYPELIST_6( T1, T2, T3, T4, T5, T6 ) scai::common::mepr::TypeList<T1,TYPELIST_5( T2, T3, T4, T5, T6 ) >
#define TYPELIST_7( T1, T2, T3, T4, T5, T6, T7 ) scai::common::mepr::TypeList<T1,TYPELIST_6( T2, T3, T4, T5, T6, T7 ) >
#define TYPELIST_8( T1, T2, T3, T4, T5, T6, T7, T8 ) scai::common::mepr::TypeList<T1,TYPELIST_7( T2, T3, T4, T5, T6, T7, T8 ) >

#define __TYPELIST( NR, ... ) TYPELIST_##NR( __VA_ARGS__ )
#define _TYPELIST( NR, ... ) __TYPELIST( NR, __VA_ARGS__ )
#define TYPELIST( NR, ... ) _TYPELIST( NR, __VA_ARGS__ )
