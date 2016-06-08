/**
 * @file lama/mepr/MatrixStorageWrapper.hpp
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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
