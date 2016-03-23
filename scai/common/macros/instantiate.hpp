/**
 * @file instantiate.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief
 * @author eschricker
 * @date 15.03.2016
 * @since 2.0.0
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
