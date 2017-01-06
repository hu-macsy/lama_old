/**
 * @file lama/mepr/DenseStorageWrapper.hpp
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
 * @brief ToDo: Missing description in ./lama/mepr/DenseStorageWrapper.hpp
 * @author eschricker
 * @date 16.03.2016
 */

#pragma once

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai
{

namespace lama
{

namespace mepr
{

template<typename ValueType, typename TList> struct DenseStorageWrapper;

template<typename ValueType>
struct DenseStorageWrapper<ValueType, common::mepr::NullType>
{
    static bool assignImpl( DenseStorage<ValueType>&, const _MatrixStorage& )
    {
        return false;
    }
};

template<typename ValueType, typename H, typename T>
struct DenseStorageWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static bool assignImpl( DenseStorage<ValueType>& obj, const _MatrixStorage& other )
    {
        if ( other.getValueType() == common::getScalarType<H>() )
        {
            const DenseStorage<H>& otherTyped = reinterpret_cast<const DenseStorage<H>& >( other );
            obj.assignDenseStorageImpl( otherTyped );
            return true;
        }
        else
        {
            return DenseStorageWrapper<ValueType, T>::assignImpl( obj, other );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
