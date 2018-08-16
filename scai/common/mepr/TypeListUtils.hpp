/**
 * @file TypeListUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Utilities for working with TypeLists
 * @author Eric Schricker
 * @date 11.04.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

namespace scai
{

namespace common
{

namespace mepr
{

/*
 * TypeListUtils
 */
template<typename TList>
struct TypeListUtils;

template<>
struct TypeListUtils<NullType>
{
    enum { size = 0 };
};

template<typename H, typename T>
struct TypeListUtils< TypeList<H, T> >
{
    enum { size = TypeListUtils<T>::size + 1 };
};

/*
 * TypeListUtils with ValueType as additional parameter
 */
template<typename ValueType, typename TList>
struct TypeListUtilsV;

template<typename ValueType>
struct TypeListUtilsV<ValueType, NullType>
{
    enum { contains = 0 };
    enum { index = -1 };
};

template<typename ValueType, typename T>
struct TypeListUtilsV<ValueType, TypeList<ValueType, T> >
{
    enum { contains = 1 };
    enum { index = TypeListUtils< T >::size };
};

template<typename ValueType, typename H, typename T>
struct TypeListUtilsV<ValueType, TypeList<H, T> >
{
    enum { contains = TypeListUtilsV<ValueType, T>::contains };
    enum { index = TypeListUtilsV<ValueType, T>::index };
};

/*
 * TypeListUtilsVLL with ValueType and two identical lists, to get a valid ValueType
 */

template<typename VT, typename TList1, typename TList2> struct TypeListUtilsVLL;

template<typename VT, typename H, typename T>
struct TypeListUtilsVLL<VT, TypeList<H, T>, NullType>
{
    typedef H ValueType;
};

template<typename VT, typename TList1, typename H, typename T>
struct TypeListUtilsVLL<VT, TList1, TypeList<H, T> >
{
    typedef typename TypeListUtilsVLL<VT, TList1, T>::ValueType ValueType;
};

template<typename VT, typename TList1, typename T>
struct TypeListUtilsVLL<VT, TList1, TypeList<VT, T> >
{
    typedef VT ValueType;
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
