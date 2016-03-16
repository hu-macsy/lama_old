/*
 * instantiate.hpp
 *
 *  Created on: Mar 15, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/config.hpp>

/*
 * Instantiate templated classes
 */

#define SCAI_COMMON_INST_CLASS_1( _class, type )        template class COMMON_DLL_IMPORTEXPORT _class<type>;
#define SCAI_COMMON_INST_CLASS_2( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_1( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_3( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_2( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_4( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_3( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_5( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_4( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_6( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_5( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_7( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_6( _class, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_8( _class, type, ... )   template class COMMON_DLL_IMPORTEXPORT _class<type>; SCAI_COMMON_INST_CLASS_7( _class, __VA_ARGS__ )

#define __SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) SCAI_COMMON_INST_CLASS_##_cnt( _class, __VA_ARGS__ )
#define _SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) __SCAI_COMMON_INST_CLASS( _class, _cnt, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS( _class, _cnt, ... ) _SCAI_COMMON_INST_CLASS( _class, _cnt, __VA_ARGS__ )

/*
 * Instantiate templated classes which expect a templated class as template argument
 */

#define SCAI_COMMON_INST_CLASS_II_1( _class, _sub, type )       template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >;
#define SCAI_COMMON_INST_CLASS_II_2( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_1( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_3( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_2( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_4( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_3( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_5( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_4( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_6( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_5( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_7( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_6( _class, _sub, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II_8( _class, _sub, type, ... )  template class COMMON_DLL_IMPORTEXPORT _class<_sub<type> >; SCAI_COMMON_INST_CLASS_II_7( _class, _sub, __VA_ARGS__ )

#define __SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) SCAI_COMMON_INST_CLASS_II_##_cnt( _class, _sub, __VA_ARGS__ )
#define _SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) __SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, __VA_ARGS__ )
#define SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, ... ) _SCAI_COMMON_INST_CLASS_II( _class, _sub, _cnt, __VA_ARGS__ )
