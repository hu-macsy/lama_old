/**
 * @file unique_name.hpp
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
 * @brief Some macro utilities for generating symbol names.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

/** Help macro to concatenate two symbols, can also deal with nested calls. */

#define SCAI_COMMON_JOIN( symbol1, symbol2 ) _SCAI_COMMON_DO_JOIN( symbol1, symbol2 )

/** Help macro to deal with nested calls of SCAI_COMMON_JOIN */

#define _SCAI_COMMON_DO_JOIN( symbol1, symbol2 ) _SCAI_COMMON_DO_JOIN2( symbol1, symbol2 )

/** Furtherhelp macro to deal with nested calls of SCAI_COMMON_JOIN */

#define _SCAI_COMMON_DO_JOIN2( symbol1, symbol2 ) symbol1##symbol2

/** @brief Creates a unique symbol name by joining the prefix, the line and the postfix.
 *
 *  \code
 *  SCAI_COMMON_UNIQUE_NAME( Interface, Registry ) -> Interface17Registry
 *  \endcode
 */

#define SCAI_COMMON_UNIQUE_NAME( prefix, postfix )                                    \
    SCAI_COMMON_JOIN( prefix , SCAI_COMMON_JOIN( __LINE__ , postfix ) )
