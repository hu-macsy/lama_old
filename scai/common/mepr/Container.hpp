#pragma once

namespace scai {

namespace common {

namespace mepr {

template<template<typename> class R, typename T1, typename T2=T1, typename T3=T2, typename T4=T3, typename T5=T4, typename T6=T5>
class Container
{
    typedef T1 head;
    typedef T6 tail;
};

template<template<typename> class R, typename T1>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1> )
{
    R<T1>::initAndReg( flag );
}

template<template<typename> class R, typename T1, typename T2>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1, T2> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, Container<R, T2>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1, T2, T3> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, Container<R, T2, T3>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1, T2, T3, T4> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, Container<R, T2, T3, T4>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1, T2, T3, T4, T5> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, Container<R, T2, T3, T4, T5>() );
}

template<template<typename> class R, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
static void instantiate( kregistry::KernelRegistry::KernelRegistryFlag flag, Container<R, T1, T2, T3, T4, T5, T6> )
{
    R<T1>::initAndReg( flag );
    instantiate( flag, Container<R, T2, T3, T4, T5, T6>() );
}

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
