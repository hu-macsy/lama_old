/**
 * @file common/mepr/TypeList.hpp
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
 * @brief TypeList structure
 * @author Eric Schricker
 * @date 03.03.2016
 */

#pragma once

namespace scai
{

namespace common
{

namespace mepr
{

/*
 * NullType, used for termination
 */
class NullType {};

template<typename H, typename T>
struct TypeList
{
    typedef H head;
    typedef T tail;
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */

/*
 * Internal used macros for TypeList creation
 *
 * TYPELIST( NR, ... ) for usage. NR is the number of types
 */

#define TYPELIST_1( T1 ) scai::common::mepr::TypeList<T1,common::mepr::NullType>
#define TYPELIST_2( T1, T2 ) scai::common::mepr::TypeList<T1,TYPELIST_1( T2 ) >
#define TYPELIST_3( T1, T2, T3 ) scai::common::mepr::TypeList<T1,TYPELIST_2( T2, T3 ) >
#define TYPELIST_4( T1, T2, T3, T4 ) scai::common::mepr::TypeList<T1,TYPELIST_3( T2, T3, T4 ) >
#define TYPELIST_5( T1, T2, T3, T4, T5 ) scai::common::mepr::TypeList<T1,TYPELIST_4( T2, T3, T4, T5 ) >
#define TYPELIST_6( T1, T2, T3, T4, T5, T6 ) scai::common::mepr::TypeList<T1,TYPELIST_5( T2, T3, T4, T5, T6 ) >
#define TYPELIST_7( T1, T2, T3, T4, T5, T6, T7 ) scai::common::mepr::TypeList<T1,TYPELIST_6( T2, T3, T4, T5, T6, T7 ) >
#define TYPELIST_8( T1, T2, T3, T4, T5, T6, T7, T8 ) scai::common::mepr::TypeList<T1,TYPELIST_7( T2, T3, T4, T5, T6, T7, T8 ) >

#define __TYPELIST( NR, ... ) TYPELIST_##NR( __VA_ARGS__ )
#define _TYPELIST( NR, ... ) __TYPELIST( NR, __VA_ARGS__ )
#define TYPELIST( NR, ... ) _TYPELIST( NR, __VA_ARGS__ )
