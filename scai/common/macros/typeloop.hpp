/*
 * loop.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/config.hpp>

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

