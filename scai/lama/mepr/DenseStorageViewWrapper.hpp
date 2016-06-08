/**
 * @file lama/mepr/DenseStorageViewWrapper.hpp
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
 * @endlicense
 *
 * @brief ToDo: Missing description in ./lama/mepr/DenseStorageViewWrapper.hpp
 * @author eschricker
 * @date 16.03.2016
 */

#pragma once

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace lama {

namespace mepr {

template<typename ValueType, typename TList> struct DenseStorageViewWrapper;

template<typename ValueType>
struct DenseStorageViewWrapper<ValueType, common::mepr::NullType>
{
    static bool assignImpl( DenseStorageView<ValueType>&, const _MatrixStorage& )
    {
        return false;
    }
};

template<typename ValueType, typename H, typename T>
struct DenseStorageViewWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static bool assignImpl( DenseStorageView<ValueType>& obj, const _MatrixStorage& other )
    {
        if( other.getValueType() == common::getScalarType<H>() )
        {
            const DenseStorageView<H>& otherTyped = reinterpret_cast<const DenseStorageView<H>& >( other );
            obj.assignDenseStorageImpl( otherTyped );
            return true;
        }
        else
        {
            return DenseStorageViewWrapper<ValueType, T>::assignImpl( obj, other );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
