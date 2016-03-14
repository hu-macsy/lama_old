/*
 * TemplateSpecifier.hpp
 *
 *  Created on: Mar 8, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

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
    static void set( FunctionType ) {}
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
