/*
 * DenseStorageWrapper.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: eschricker
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
