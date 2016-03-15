/**
 * @file lama/mepr/MatrixStorageWrapper.hpp
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
 * @brief Wrapper for templated function calls in MatrixStorage
 * @author Eric Schricker
 * @date 14.03.2016
 */

#pragma once

#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/mepr/TypeList.hpp>
#include <scai/hmemo/HArray.hpp>

namespace scai {

namespace lama {

namespace mepr {

/*
 * Forward declaration
 */
template<typename ValueType, typename TList>
struct MatrixStorageWrapper;

/*
 * Termination
 */
template<typename ValueType>
struct MatrixStorageWrapper<ValueType, common::mepr::NullType>
{
    static void setDenseData( MatrixStorage<ValueType>*, const IndexType, const IndexType, const hmemo::_HArray&, const ValueType )
    { }
};

/*
 * Step n
 */
template<typename ValueType, typename H, typename T>
struct MatrixStorageWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static void setDenseData(
        MatrixStorage<ValueType>* obj,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::_HArray& values,
        const ValueType epsilon )
    {
        if( values.getValueType() == common::getScalarType<H>() )
        {
            hmemo::_HArray& mValues = const_cast<hmemo::_HArray&>( values );
            hmemo::HArray<H>& typedValues = reinterpret_cast<hmemo::HArray<H>& >( mValues );
            const DenseStorageView<H> denseStorage( typedValues, numRows, numColumns );
            H tmpEpsilon = static_cast<H>( epsilon );
            denseStorage.swapEpsilon( tmpEpsilon );
            obj->assign( denseStorage );
        }
        else
        {
            MatrixStorageWrapper<ValueType, T>::setDenseData( obj, numRows, numColumns, values, epsilon );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
