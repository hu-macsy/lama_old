/**
 * @file common/mepr/TemplateSpecifier.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
