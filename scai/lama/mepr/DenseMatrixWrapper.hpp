/*
 * DenseMatrixWrapper.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace lama {

namespace mepr {

template<typename ValueType, typename TList>
struct DenseMatrixWrapper;

template<typename ValueType>
struct DenseMatrixWrapper<ValueType, common::mepr::NullType>
{
    static void assignDenseImpl( DenseMatrix<ValueType>&, const Matrix& other)
    {
        COMMON_THROWEXCEPTION( "type dense matrix not supported --> " << other )
    }
    static void assignSparseImpl( DenseMatrix<ValueType>&, const Matrix& other )
    {
        COMMON_THROWEXCEPTION( "type of sparse matrix not supported --> " << other )
    }

    static void getDiagonalImpl( const DenseMatrix<ValueType>&, Vector& other )
    {
        COMMON_THROWEXCEPTION( "type of vector not supported --> " << other )
    }
};

template<typename ValueType, typename H, typename T>
struct DenseMatrixWrapper<ValueType, common::mepr::TypeList<H,T> >
{
    static void assignDenseImpl( DenseMatrix<ValueType>& obj, const Matrix& other)
    {
        if( other.getValueType() == common::getScalarType<H>() )
        {
            obj.copyDenseMatrix( reinterpret_cast<const DenseMatrix<H>& >( other ) );
        }
        else
        {
            DenseMatrixWrapper<ValueType, T>::assignDenseImpl( obj, other );
        }
    }

    static void assignSparseImpl( DenseMatrix<ValueType>& obj, const Matrix& other )
    {
        if( other.getValueType() == common::getScalarType<H>() )
        {
            const SparseMatrix<H>& sparse = reinterpret_cast<const SparseMatrix<H>& >( other );
            const CSRSparseMatrix<ValueType> tmp = sparse;
            obj.assignSparse( tmp );
        }
        else
        {
            DenseMatrixWrapper<ValueType, T>::assignSparseImpl( obj, other );
        }
    }

    static void getDiagonalImpl( const DenseMatrix<ValueType>& obj, Vector& diagonal )
    {
        if( diagonal.getValueType() == common::getScalarType<H>() )
        {
            DenseVector<H>& dense = reinterpret_cast<DenseVector<H>& >( diagonal );
            obj.getDiagonalImpl( dense );
        }
        else
        {
            DenseMatrixWrapper<ValueType, T>::getDiagonalImpl( obj, diagonal );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
