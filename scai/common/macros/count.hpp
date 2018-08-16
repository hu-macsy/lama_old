/**
 * @file count.hpp
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
 * @brief Definition of macro that counts variadic macro arguments
 * @author Thomas Brandes
 * @date 20.05.2016
 */

#pragma once

/** This macro counts the number of arguments for a variadic list
 *
 *  \code
 *     SCAI_COMMON_COUT_NARG( A, B, C, D, E ) -> 5
 *     SCAI_COMMON_COUT_NARG( ) -> 0
 *  \endcode
 */

#define SCAI_COMMON_COUNT_NARG( ... ) COUNT_NARG_( __VA_ARGS__ )
#define COUNT_NARG_( ... ) COUNT_ARGS_( , ##__VA_ARGS__, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define COUNT_ARGS_( z, a, b, c, d, e, f, g, h, i, j, k, l, cnt, ...) cnt

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

