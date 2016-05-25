/**
 * @file kregistry/mepr/Registrator.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Registrator structure, for registering templated functions in KernelRegistry
 * @author Eric Schricker
 * @date 03.03.2016
 */

#pragma once

// library
#include <scai/kregistry/KernelRegistry.hpp>


/**
 * Declare Registrator inside a class
 *
 * Functions without template-parameter
 * - SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
 *
 * Functions with one template-parameter
 * - SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
 *
 * Functions with two template-parameter
 * - SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )
 *
 */

#define SCAI_KREGISTRY_DECL_REGISTRATOR( name, ... )                                                     \
    __VA_ARGS__                                                                                   \
    struct name                                                                                   \
    {                                                                                             \
        static void initAndReg( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag ); \
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
    static void call( const KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename> class R, typename H, typename T>
struct RegistratorV< R, common::mepr::TypeList<H,T> >
{
    static void call( const KernelRegistry::KernelRegistryFlag flag )
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
    static void call( const KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename ValueType>
struct _RegistratorVO<R, ValueType, common::mepr::NullType>
{
    static void call( const KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename ValueType, typename H, typename T>
struct _RegistratorVO< R, ValueType, common::mepr::TypeList<H,T> >
{
    static void call( const KernelRegistry::KernelRegistryFlag flag )
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
    static void call( const KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename TList >
struct RegistratorVO<R, common::mepr::NullType, TList >
{
    static void call( const KernelRegistry::KernelRegistryFlag ){}
};

template<template<typename,typename> class R, typename H1, typename T1, typename H2, typename T2>
struct RegistratorVO< R, common::mepr::TypeList<H1,T1>, common::mepr::TypeList<H2, T2> >
{
    static void call( const KernelRegistry::KernelRegistryFlag flag )
    {
        _RegistratorVO<R, H1, common::mepr::TypeList<H2,T2> >::call( flag );

        RegistratorVO<R, T1, common::mepr::TypeList<H2,T2> >::call( flag );
    }
};

} /* end namespace mepr */

} /* end namespace kregistry */

} /* end namespace scai */
