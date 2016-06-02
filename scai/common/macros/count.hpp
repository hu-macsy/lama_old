/**
 * @file count.hpp
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
 * @brief Definition of macro that counts variadic macro arguments
 * @author Thomas Brandes
 * @date 20.05.2016
 */

#pragma once

#define SCAI_COMMON_COUNT_ARG_N( \
          _1,  _2,  _3,  _4,  _5,  _6,  _7,  _8,  _9, _10, \
         _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, \
         _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, \
         _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, \
         _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, \
         _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, \
         _61, _62, _63, N, ...) N

#define SCAI_COMMON_COUNT_RSEQ_N()      63, 62, 61, 60, \
         59, 58, 57, 56, 55, 54, 53, 52, 51, 50, \
         49, 48, 47, 46, 45, 44, 43, 42, 41, 40, \
         39, 38, 37, 36, 35, 34, 33, 32, 31, 30, \
         29, 28, 27, 26, 25, 24, 23, 22, 21, 20, \
         19, 18, 17, 16, 15, 14, 13, 12, 11, 10, \
          9,  8,  7,  6,  5,  4,  3,  2,  1,  0

#define SCAI_COMMON_COUNT_NARG_( ... ) SCAI_COMMON_COUNT_ARG_N( __VA_ARGS__ )

/** This macro counts the number of arguments for a variadic list
 *
 *  \code
 *     SCAI_COMMON_COUT_NARG( A, B, C, D, E ) -> 5
 *     SCAI_COMMON_COUT_NARG( ) -> 0 
 *  \endcode
 */

#define SCAI_COMMON_COUNT_NARG(...) SCAI_COMMON_COUNT_NARG_( __VA_ARGS__, SCAI_COMMON_COUNT_RSEQ_N() )

/** Macro gives the first argument of a variadic argument list. Uses help macro to force evaluation of macro before
 *
 *  \code
 *     SCAI_COMMON_FIRST_ARG( A, B, C, D, E ) -> A
 *  \endcode
 */

#define SCAI_COMMON_FIRST_ARG( x, ... ) _SCAI_COMMON_FIRST_ARG( x, __VA_ARGS__ )

/** Help macro needed as x must be replaced first */

#define _SCAI_COMMON_FIRST_ARG( x, ... ) x

/** Macro gives the tail of a variadic argument list, first removed
 *
 *  \code
 *     SCAI_COMMON_TAIL( A, B, C, D, E ) -> A
 *  \endcode
 */

#define SCAI_COMMON_TAIL( x, ... ) _SCAI_COMMON_TAIL( x, __VA_ARGS__ )

/** Help macro needed as x must be replaced first */

#define _SCAI_COMMON_TAIL( x, ... ) __VA_ARGS__

