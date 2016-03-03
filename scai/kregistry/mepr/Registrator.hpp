#pragma once

// library
#include <scai/kregistry/KernelRegistry.hpp>

// internal scai libraries
#include <scai/common/mepr/Container.hpp>

/*
 * Declare Registrator inside a class
 *
 * Functions without template-parameter
 * - SCAI_DECLARE_REGISTRATOR( Registrator )
 *
 * Functions with one template-parameter
 * - SCAI_DECLARE_REGISTRATOR( RegistratorV, template<typename ValueType> )
 *
 * Functions with two template-parameter
 * - SCAI_DECLARE_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )
 *
 */

#define SCAI_DECLARE_REGISTRATOR( name, ... )                                               \
    __VA_ARGS__                                                                             \
    struct name                                                                             \
    {                                                                                       \
        static void initAndReg( const kregistry::KernelRegistry::KernelRegistryFlag flag ); \
    };

namespace scai {

namespace kregistry {

namespace mepr {

/*
 * Instantiate the Registrator for one template-parameter with different kind of ValueTypes
 */

template<template<typename> class R, typename TList> struct RegistratorV;

template<template<typename> class R> struct RegistratorV<R,common::mepr::NullType>
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename> class R, typename H, typename T>
struct RegistratorV< R, common::mepr::TypeList<H,T> >
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag flag )
    {
        R<H>::initAndReg( flag );
        RegistratorV<R,T>::call( flag );
    }
};

/*
 * Instantiate the Registrator for two template-parameter with different kind of ValueTypes
 * Every combination will be instantiated
 *
 * _RegistratorVO is just internal used
 */

/*
 * _RegistratorVO
 */

template<template<typename,typename> class R, typename ValueType, typename TList> struct _RegistratorVO;

template<template<typename,typename> class R, typename TList>
struct _RegistratorVO<R, common::mepr::NullType, TList>
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename ValueType>
struct _RegistratorVO<R, ValueType, common::mepr::NullType>
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename ValueType, typename H, typename T>
struct _RegistratorVO< R, ValueType, common::mepr::TypeList<H,T> >
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag flag )
    {
        R<ValueType, H>::initAndReg( flag );
        _RegistratorVO<R, ValueType, T>::call( flag );
    }
};

/*
 * RegistratorVO
 */

template<template<typename,typename> class R, typename TList1, typename TList2> struct RegistratorVO;

template<template<typename,typename> class R >
struct RegistratorVO<R, common::mepr::NullType, common::mepr::NullType >
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename TList >
struct RegistratorVO<R, common::mepr::NullType, TList >
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename H1, typename T1, typename H2, typename T2>
struct RegistratorVO< R, common::mepr::TypeList<H1,T1>, common::mepr::TypeList<H2, T2> >
{
    static void call( const kregistry::KernelRegistry::KernelRegistryFlag flag )
    {
        _RegistratorVO<R, H1, common::mepr::TypeList<H2,T2> >::call( flag );

        RegistratorVO<R, T1, common::mepr::TypeList<H2,T2> >::call( flag );
    }
};

} /* end namespace mepr */

} /* end namespace kregistry */

} /* end namespace scai */
