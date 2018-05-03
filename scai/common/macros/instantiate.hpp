/**
 * @file macros/instantiate.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Macro for template instantiation                      
 * @author eschricker
 * @date 15.03.2016
 */

#pragma once

#include <scai/common/config.hpp>

/*
 * Instantiate templated classes
 */

#define SCAI_COMMON_INST_CLASS_1( _class, type )        template class _class<type>;
#define SCAI_COMMON_INST_CLASS_2( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_1( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_3( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_2( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_4( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_3( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_5( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_4( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_6( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_5( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_7( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_6( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_8( _class, type, ... )   template class _class<type>; SCAI_COMMON_INST_CLASS_7( _class, __VA_ARGS__ )

#define ___SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) SCAI_COMMON_INST_CLASS_##_cnt( _class, __VA_ARGS__ )
#define __SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) ___SCAI_COMMON_INST_CLASS( _class, _cnt, __VA_ARGS__ )
#define _SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) __SCAI_COMMON_INST_CLASS( _class, _cnt, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS( _class, ... ) _SCAI_COMMON_INST_CLASS( _class, SCAI_COMMON_COUNT_NARG( __VA_ARGS__ ), __VA_ARGS__ )

/*
 * Instantiate templated classes which expect a templated class as template argument
 */

#define SCAI_COMMON_INST_CLASS_II_1( _class, _sub, type )       template class _class<_sub<type> >;
#define SCAI_COMMON_INST_CLASS_II_2( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_1( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_3( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_2( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_4( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_3( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_5( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_4( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_6( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_5( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_7( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_6( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_8( _class, _sub, type, ... )  template class _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_7( _class, _sub, __VA_ARGS__ )

#define ___SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) SCAI_COMMON_INST_CLASS_II_##_cnt( _class, _sub, __VA_ARGS__ )
#define __SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) ___SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, __VA_ARGS__ )
#define _SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) __SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II( _class, _sub, ... ) _SCAI_COMMON_INST_CLASS_II( _class, _sub, SCAI_COMMON_COUNT_NARG( __VA_ARGS__), __VA_ARGS__ )
