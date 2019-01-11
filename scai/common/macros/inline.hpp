/**
 * @file inline.hpp
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
 * @brief Macro for static inline routines handled differently by compilers.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

/** Macros for inlining used on different OS */

#ifdef WIN32
#define LAMA_STATIC_INLINE_FUNCTION_PREFIX static __inline
#define LAMA_INLINE_FUNCTION_PREFIX __inline
#else
#define LAMA_STATIC_INLINE_FUNCTION_PREFIX static inline
#define LAMA_INLINE_FUNCTION_PREFIX static inline
#endif
