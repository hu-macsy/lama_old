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
