/**
 * @file macros/loop.hpp
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
 * @brief Definition of macro that applies another macro to a variadic argument list
 * @author eschricker
 * @date 16.03.2016
 */

#pragma once

#include <scai/common/config.hpp>
#include <scai/common/macros/count.hpp>

/*
 * Level 1
 */

#define SCAI_COMMON_LOOP_0( _macro, type )
#define SCAI_COMMON_LOOP_1( _macro, type ) _macro( type )
#define SCAI_COMMON_LOOP_2( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_1( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_3( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_2( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_4( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_3( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_5( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_4( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_6( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_5( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_7( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_6( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_8( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_7( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_9( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_8( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_10( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_9( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_11( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_10( _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_12( _macro, type, ... ) _macro( type ) SCAI_COMMON_LOOP_11( _macro, __VA_ARGS__ )

#define ___SCAI_COMMON_LOOP( _cnt, _macro, ... ) SCAI_COMMON_LOOP_##_cnt( _macro, __VA_ARGS__ )
#define __SCAI_COMMON_LOOP( _cnt, _macro, ... ) ___SCAI_COMMON_LOOP( _cnt, _macro, __VA_ARGS__ )
#define _SCAI_COMMON_LOOP( _cnt, _macro, ... ) __SCAI_COMMON_LOOP( _cnt, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP( _macro, ... ) __SCAI_COMMON_LOOP( SCAI_COMMON_COUNT_NARG( __VA_ARGS__ ), _macro, __VA_ARGS__ )

/*
 * Level 2
 */

#define SCAI_COMMON_LOOP_LVL2_1( arg1, _macro, type ) _macro( arg1, type )
#define SCAI_COMMON_LOOP_LVL2_2( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_1( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_3( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_2( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_4( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_3( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_5( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_4( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_6( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_5( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_7( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_6( arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2_8( arg1, _macro, type, ... ) _macro( arg1, type ) SCAI_COMMON_LOOP_LVL2_7( arg1, _macro, __VA_ARGS__ )

#define __SCAI_COMMON_LOOP_LVL2( _cnt, arg1, _macro, ... ) SCAI_COMMON_LOOP_LVL2_##_cnt( arg1, _macro, __VA_ARGS__ )
#define _SCAI_COMMON_LOOP_LVL2( _cnt, arg1, _macro, ... ) __SCAI_COMMON_LOOP_LVL2( _cnt, arg1, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL2( arg1, _macro, ... ) _SCAI_COMMON_LOOP_LVL2( SCAI_COMMON_COUNT_NARG( __VA_ARGS__ ), arg1, _macro, __VA_ARGS__ )

/*
 * Level 3
 */

#define SCAI_COMMON_LOOP_LVL3_1( arg1, arg2, _macro, type ) _macro( arg1, arg2, type )
#define SCAI_COMMON_LOOP_LVL3_2( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_1( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_3( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_2( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_4( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_3( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_5( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_4( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_6( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_5( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_7( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_6( arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3_8( arg1, arg2, _macro, type, ... ) _macro( arg1, arg2, type ) SCAI_COMMON_LOOP_LVL3_7( arg1, arg2, _macro, __VA_ARGS__ )

#define __SCAI_COMMON_LOOP_LVL3( _cnt, arg1, arg2, _macro, ... ) SCAI_COMMON_LOOP_LVL3_##_cnt( arg1, arg2, _macro, __VA_ARGS__ )
#define _SCAI_COMMON_LOOP_LVL3( _cnt, arg1, arg2, _macro, ... ) __SCAI_COMMON_LOOP_LVL3( _cnt, arg1, arg2, _macro, __VA_ARGS__ )
#define SCAI_COMMON_LOOP_LVL3( arg1, arg2, _macro, ... ) _SCAI_COMMON_LOOP_LVL3( SCAI_COMMON_COUNT_NARG( __VA_ARGS__ ), arg1, arg2, _macro, __VA_ARGS__ )

