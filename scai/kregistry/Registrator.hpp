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
        static void initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag );       \
    };

namespace scai {

namespace kregistry {

/*
 * Instantiate the Registrator for one template-parameter with different kind of ValueTypes
 * Currently up to 12 supported
 */

template<template<typename> class R, typename T1>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1> )
{
    R<T1>::initAndReg( flag );
}

template<template<typename> class R, typename T1, typename T2>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7, T8> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7, T8>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7, T8, T9> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7, T8, T9>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7, T8, T9, T10>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10, typename T11>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerV<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, common::mepr::ContainerV<R, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12>() );
}

/*
 * Instantiate the Registrator for two template-parameter with different kind of ValueTypes
 * Currently up to 12 supported
 * Every combination will be instantiated
 */

template<template<typename, typename> class R, typename T1, typename T2>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2> )
{
    R<T1, T1>::initAndReg( flag );
    R<T1, T2>::initAndReg( flag );
    R<T2, T1>::initAndReg( flag );
    R<T2, T2>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2>() );

    R<T1, T3>::initAndReg( flag );
    R<T2, T3>::initAndReg( flag );

    R<T3, T1>::initAndReg( flag );
    R<T3, T2>::initAndReg( flag );
    R<T3, T3>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3>() );

    R<T1, T4>::initAndReg( flag );
    R<T2, T4>::initAndReg( flag );
    R<T3, T4>::initAndReg( flag );

    R<T4, T1>::initAndReg( flag );
    R<T4, T2>::initAndReg( flag );
    R<T4, T3>::initAndReg( flag );
    R<T4, T4>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4>() );

    R<T1, T5>::initAndReg( flag );
    R<T2, T5>::initAndReg( flag );
    R<T3, T5>::initAndReg( flag );
    R<T4, T5>::initAndReg( flag );

    R<T5, T1>::initAndReg( flag );
    R<T5, T2>::initAndReg( flag );
    R<T5, T3>::initAndReg( flag );
    R<T5, T4>::initAndReg( flag );
    R<T5, T5>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5>() );

    R<T1, T6>::initAndReg( flag );
    R<T2, T6>::initAndReg( flag );
    R<T3, T6>::initAndReg( flag );
    R<T4, T6>::initAndReg( flag );
    R<T5, T6>::initAndReg( flag );

    R<T6, T1>::initAndReg( flag );
    R<T6, T2>::initAndReg( flag );
    R<T6, T3>::initAndReg( flag );
    R<T6, T4>::initAndReg( flag );
    R<T6, T5>::initAndReg( flag );
    R<T6, T6>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5,
        typename T6, typename T7>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6>() );

    R<T1, T7>::initAndReg( flag );
    R<T2, T7>::initAndReg( flag );
    R<T3, T7>::initAndReg( flag );
    R<T4, T7>::initAndReg( flag );
    R<T5, T7>::initAndReg( flag );
    R<T6, T7>::initAndReg( flag );

    R<T7, T1>::initAndReg( flag );
    R<T7, T2>::initAndReg( flag );
    R<T7, T3>::initAndReg( flag );
    R<T7, T4>::initAndReg( flag );
    R<T7, T5>::initAndReg( flag );
    R<T7, T6>::initAndReg( flag );
    R<T7, T7>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5,
        typename T6, typename T7, typename T8>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7>() );

    R<T1, T8>::initAndReg( flag );
    R<T2, T8>::initAndReg( flag );
    R<T3, T8>::initAndReg( flag );
    R<T4, T8>::initAndReg( flag );
    R<T5, T8>::initAndReg( flag );
    R<T6, T8>::initAndReg( flag );
    R<T7, T8>::initAndReg( flag );

    R<T8, T1>::initAndReg( flag );
    R<T8, T2>::initAndReg( flag );
    R<T8, T3>::initAndReg( flag );
    R<T8, T4>::initAndReg( flag );
    R<T8, T5>::initAndReg( flag );
    R<T8, T6>::initAndReg( flag );
    R<T8, T7>::initAndReg( flag );
    R<T8, T8>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8>() );

    R<T1, T9>::initAndReg( flag );
    R<T2, T9>::initAndReg( flag );
    R<T3, T9>::initAndReg( flag );
    R<T4, T9>::initAndReg( flag );
    R<T5, T9>::initAndReg( flag );
    R<T6, T9>::initAndReg( flag );
    R<T7, T9>::initAndReg( flag );
    R<T8, T9>::initAndReg( flag );

    R<T9, T1>::initAndReg( flag );
    R<T9, T2>::initAndReg( flag );
    R<T9, T3>::initAndReg( flag );
    R<T9, T4>::initAndReg( flag );
    R<T9, T5>::initAndReg( flag );
    R<T9, T6>::initAndReg( flag );
    R<T9, T7>::initAndReg( flag );
    R<T9, T8>::initAndReg( flag );
    R<T9, T9>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9>() );

    R<T1, T10>::initAndReg( flag );
    R<T2, T10>::initAndReg( flag );
    R<T3, T10>::initAndReg( flag );
    R<T4, T10>::initAndReg( flag );
    R<T5, T10>::initAndReg( flag );
    R<T6, T10>::initAndReg( flag );
    R<T7, T10>::initAndReg( flag );
    R<T8, T10>::initAndReg( flag );
    R<T9, T10>::initAndReg( flag );

    R<T10, T1>::initAndReg( flag );
    R<T10, T2>::initAndReg( flag );
    R<T10, T3>::initAndReg( flag );
    R<T10, T4>::initAndReg( flag );
    R<T10, T5>::initAndReg( flag );
    R<T10, T6>::initAndReg( flag );
    R<T10, T7>::initAndReg( flag );
    R<T10, T8>::initAndReg( flag );
    R<T10, T9>::initAndReg( flag );
    R<T10, T10>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10, typename T11>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>() );

    R<T1, T11>::initAndReg( flag );
    R<T2, T11>::initAndReg( flag );
    R<T3, T11>::initAndReg( flag );
    R<T4, T11>::initAndReg( flag );
    R<T5, T11>::initAndReg( flag );
    R<T6, T11>::initAndReg( flag );
    R<T7, T11>::initAndReg( flag );
    R<T8, T11>::initAndReg( flag );
    R<T9, T11>::initAndReg( flag );
    R<T10, T11>::initAndReg( flag );

    R<T11, T1>::initAndReg( flag );
    R<T11, T2>::initAndReg( flag );
    R<T11, T3>::initAndReg( flag );
    R<T11, T4>::initAndReg( flag );
    R<T11, T5>::initAndReg( flag );
    R<T11, T6>::initAndReg( flag );
    R<T11, T7>::initAndReg( flag );
    R<T11, T8>::initAndReg( flag );
    R<T11, T9>::initAndReg( flag );
    R<T11, T10>::initAndReg( flag );
    R<T11, T11>::initAndReg( flag );
}

template<template<typename, typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
    typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
static void instantiate(
    kregistry::KernelRegistry::KernelRegistryFlag flag,
    common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> )
{
    instantiate( flag, common::mepr::ContainerVO<R, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>() );

    R<T1, T12>::initAndReg( flag );
    R<T2, T12>::initAndReg( flag );
    R<T3, T12>::initAndReg( flag );
    R<T4, T12>::initAndReg( flag );
    R<T5, T12>::initAndReg( flag );
    R<T6, T12>::initAndReg( flag );
    R<T7, T12>::initAndReg( flag );
    R<T8, T12>::initAndReg( flag );
    R<T9, T12>::initAndReg( flag );
    R<T10, T12>::initAndReg( flag );
    R<T11, T12>::initAndReg( flag );

    R<T12, T1>::initAndReg( flag );
    R<T12, T2>::initAndReg( flag );
    R<T12, T3>::initAndReg( flag );
    R<T12, T4>::initAndReg( flag );
    R<T12, T5>::initAndReg( flag );
    R<T12, T6>::initAndReg( flag );
    R<T12, T7>::initAndReg( flag );
    R<T12, T8>::initAndReg( flag );
    R<T12, T9>::initAndReg( flag );
    R<T12, T10>::initAndReg( flag );
    R<T12, T11>::initAndReg( flag );
    R<T12, T12>::initAndReg( flag );
}

} /* end namespace kregistry */

} /* end namespace scai */
