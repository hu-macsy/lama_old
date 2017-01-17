/**
 * @file lama/mepr/DenseMatrixWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief ToDo: Missing description in ./lama/mepr/DenseMatrixWrapper.hpp
 * @author eschricker
 * @date 16.03.2016
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai
{

namespace lama
{

namespace mepr
{

template<typename ValueType, typename TList>
struct DenseMatrixWrapper;

template<typename ValueType>
struct DenseMatrixWrapper<ValueType, common::mepr::NullType>
{
    static void assignDenseImpl( DenseMatrix<ValueType>&, const Matrix& other )
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
struct DenseMatrixWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static void assignDenseImpl( DenseMatrix<ValueType>& obj, const Matrix& other )
    {
        if ( other.getValueType() == common::getScalarType<H>() )
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
        if ( other.getValueType() == common::getScalarType<H>() )
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
        if ( diagonal.getValueType() == common::getScalarType<H>() )
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
