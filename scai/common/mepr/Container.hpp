#pragma once

#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace common {

namespace mepr {

template<typename T1 = NullType, typename T2=T1, typename T3=T2, typename T4=T3, typename T5=T4, typename T6=T5,
    typename T7=T6, typename T8=T7, typename T9=T8, typename T10=T9, typename T11=T10, typename T12=T11>
class Container
{
};


/*
 * Container for classes with one template-parameter
 */

template<template<typename> class R,
    typename T1, typename T2=T1, typename T3=T2, typename T4=T3, typename T5=T4, typename T6=T5,
    typename T7=T6, typename T8=T7, typename T9=T8, typename T10=T9, typename T11=T10, typename T12=T11>
class ContainerV
{
};

/*
 * Container for classes with two template-parameters
 */

template<template<typename R1, typename R2=R1> class R,
    typename T1, typename T2=T1, typename T3=T2, typename T4=T3, typename T5=T4, typename T6=T5,
    typename T7=T6, typename T8=T7, typename T9=T8, typename T10=T9, typename T11=T10, typename T12=T11>
class ContainerVO
{
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
