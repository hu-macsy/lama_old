/**
 * @file macros/typeloop.hpp
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
 * @brief ToDo: Missing description in ./macros/typeloop.hpp
 * @author eschricker
 * @date 16.03.2016
 */

#pragma once

#include <scai/common/config.hpp>

#define __SCAI_COMMON_FIRST_ENTRY( x, ... ) x
#define _SCAI_COMMON_FIRST_ENTRY( x, ... ) __SCAI_COMMON_FIRST_ENTRY( x, __VA_ARGS__ )
#define SCAI_COMMON_FIRST_ENTRY( x, ... ) _SCAI_COMMON_FIRST_ENTRY( x, __VA_ARGS__ )

/*
 * Level 1
 */

#define SCAI_COMMON_TYPELOOP_1( _macro, type ) _macro( type )
#define SCAI_COMMON_TYPELOOP_2( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_1( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_3( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_2( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_4( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_3( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_5( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_4( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_6( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_5( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_7( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_6( _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_8( _macro, type, ... ) _macro( type ) SCAI_COMMON_TYPELOOP_7( _macro, __VA_ARGS__ )

#define __SCAI_COMMON_TYPELOOP( _cnt, _macro, ... ) SCAI_COMMON_TYPELOOP_##_cnt( _macro, __VA_ARGS__ )
#define _SCAI_COMMON_TYPELOOP( _cnt, _macro, ... ) __SCAI_COMMON_TYPELOOP( _cnt, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP( _cnt, _macro, ... ) _SCAI_COMMON_TYPELOOP( _cnt, _macro, __VA_ARGS__ )

/*
 * Level 2
 */

#define SCAI_COMMON_TYPELOOP_LVL2_1( ValueType, _macro, type ) _macro( ValueType, type )
#define SCAI_COMMON_TYPELOOP_LVL2_2( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_1( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_3( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_2( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_4( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_3( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_5( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_4( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_6( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_5( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_7( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_6( ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2_8( ValueType, _macro, type, ... ) _macro( ValueType, type ) SCAI_COMMON_TYPELOOP_LVL2_7( ValueType, _macro, __VA_ARGS__ )

#define __SCAI_COMMON_TYPELOOP_LVL2( _cnt, ValueType, _macro, ... ) SCAI_COMMON_TYPELOOP_LVL2_##_cnt( ValueType, _macro, __VA_ARGS__ )
#define _SCAI_COMMON_TYPELOOP_LVL2( _cnt, ValueType, _macro, ... ) __SCAI_COMMON_TYPELOOP_LVL2( _cnt, ValueType, _macro, __VA_ARGS__ )
#define SCAI_COMMON_TYPELOOP_LVL2( _cnt, ValueType, _macro, ... ) _SCAI_COMMON_TYPELOOP_LVL2( _cnt, ValueType, _macro, __VA_ARGS__ )

