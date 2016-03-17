/**
 * @file common/mepr/TemplateSpecifier.hpp
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
 * @brief Structure for specifing templated functions for a given TypeList
 * @author Eric Schricker
 * @date 08.03.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

#include <list>

#define SCAI_DECLARE_TEMPLATESPECIFIER( name, ... )                                         \
    __VA_ARGS__                                                                             \
    struct name                                                                             \
    {                                                                                       \
        static void specify(  );                                                            \
    };

namespace scai {

namespace common {

namespace mepr {

class TemplateSpecifier
{
public:
    template<typename FunctionType>
    static void set( FunctionType x)
    {
        theLastOne = ( void(*)() ) x;
    }

    static void (*theLastOne)();
};

/*
 * One template parameter
 *
 * TemplateSpecifierV
 */

template<template<typename> class R, typename TList> struct TemplateSpecifierV;

template<template<typename> class R> struct TemplateSpecifierV<R,common::mepr::NullType>
{
    static void call( ){}
};

template<template<typename> class R, typename H, typename T>
struct TemplateSpecifierV< R, common::mepr::TypeList<H,T> >
{
    static void call(  )
    {
        R<H>::specify( );
        TemplateSpecifierV<R,T>::call( );
    }
};

/*
 * Two template parameters
 *
 * _TemplateSpecifierVO
 *   - just used internally
 */

template<template<typename,typename> class R, typename ValueType, typename TList> struct _TemplateSpecifierVO;

template<template<typename,typename> class R, typename TList>
struct _TemplateSpecifierVO<R, common::mepr::NullType, TList>
{
    static void call( ){}
};

template<template<typename,typename> class R, typename ValueType>
struct _TemplateSpecifierVO<R, ValueType, common::mepr::NullType>
{
    static void call( ){}
};

template<template<typename,typename> class R, typename ValueType, typename H, typename T>
struct _TemplateSpecifierVO< R, ValueType, common::mepr::TypeList<H,T> >
{
    static void call( )
    {
        R<ValueType, H>::specify( );
        _TemplateSpecifierVO<R, ValueType, T>::call( );
    }
};

/*
 * TemplateSpecifierVO
 *   -  used for template specification
 */

template<template<typename,typename> class R, typename TList1, typename TList2> struct TemplateSpecifierVO;

template<template<typename,typename> class R >
struct TemplateSpecifierVO<R, common::mepr::NullType, common::mepr::NullType >
{
    static void call( ){}
};

template<template<typename,typename> class R, typename TList >
struct TemplateSpecifierVO<R, common::mepr::NullType, TList >
{
    static void call( ){}
};

template<template<typename,typename> class R, typename H1, typename T1, typename H2, typename T2>
struct TemplateSpecifierVO< R, common::mepr::TypeList<H1,T1>, common::mepr::TypeList<H2, T2> >
{
    static void call( )
    {
        _TemplateSpecifierVO<R, H1, common::mepr::TypeList<H2,T2> >::call( );

        TemplateSpecifierVO<R, T1, common::mepr::TypeList<H2,T2> >::call( );
    }
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
