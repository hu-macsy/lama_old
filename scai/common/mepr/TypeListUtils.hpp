/**
 * @file TypeListUtils.hpp
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
 * @brief Utilities for working with TypeLists
 * @author Eric Schricker
 * @date 11.04.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace common
{

namespace mepr
{

/*
 * TemplateListUtils
 */
template<typename TList>
struct TypeListUtils;

template<>
struct TypeListUtils<NullType>
{
    enum{ size = 0 };
};

template<typename H, typename T>
struct TypeListUtils< TypeList<H, T> >
{
    enum{ size = TypeListUtils<T>::size + 1 };
};

/*
 * TypeListUtils with ValueType as additional parameter
 */
template<typename ValueType, typename TList>
struct TypeListUtilsV;

template<typename ValueType>
struct TypeListUtilsV<ValueType, NullType>
{
    enum{ contains = 0 };
    enum{ index = -1 };
};

template<typename ValueType, typename T>
struct TypeListUtilsV<ValueType, TypeList<ValueType,T> >
{
    enum{ contains = 1 };
    enum{ index = TypeListUtils< T >::size };
};

template<typename ValueType, typename H, typename T>
struct TypeListUtilsV<ValueType, TypeList<H,T> >
{
    enum{ contains = TypeListUtilsV<ValueType, T>::contains };
    enum{ index = TypeListUtilsV<ValueType, T>::index };
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
