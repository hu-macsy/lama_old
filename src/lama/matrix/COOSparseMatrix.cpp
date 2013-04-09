/**
 * @file COOSparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of methods and constructors for template class COOSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 * $Id$
 */

// hpp
#include <lama/matrix/COOSparseMatrix.hpp>

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, COOSparseMatrix<T>::logger, "Matrix.SparseMatrix.COOSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    LAMA_LOG_INFO( logger, "COOSpareMatrix()" )
}

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    LAMA_LOG_INFO( logger, "COOSpareMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::~COOSparseMatrix()
{
    LAMA_LOG_INFO( logger, "~COOSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=( const COOSparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, "COOSparseMatrix = COOSparseMatrix : " << matrix )
    assign( matrix );
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=( const Matrix& matrix )
{
    LAMA_LOG_INFO( logger, " = Matrix : " << matrix )
    this->assign( matrix ); // matrix does not depend on template parameter, so this-> is needed.
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=( const Expression<Matrix,Matrix,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = A * B " )
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=( const Expression<Scalar,Matrix,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = alpha * A " )
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=(
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = alpha * A * B" )
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=(
    const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    LAMA_LOG_INFO( logger, "COOSparseMatrix = alpha * A * B" )
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=(
    const Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    LAMA_LOG_INFO( logger, "COOSparseMatrix = alpha * A + beta * B" )
    Matrix::operator=( exp );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename COOSparseMatrix<ValueType>::StorageType&
COOSparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "COOSparseMatrix: local storage is no more COO: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename COOSparseMatrix<ValueType>::StorageType&
COOSparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "COOSparseMatrix: local storage is no more COO: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename COOSparseMatrix<ValueType>::StorageType&
COOSparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );

    LAMA_ASSERT_ERROR( halo, "COOSparseMatrix: halo storage is no more COO: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
{
    // make sure that local storage fits into this sparse matrix

    LAMA_ASSERT_EQUAL_ERROR( localStorage.getNumRows(), mLocalData->getNumRows() )
    LAMA_ASSERT_EQUAL_ERROR( localStorage.getNumColumns(), mLocalData->getNumColumns() )

    // make sure that local matrix storage has the correct format / value type

    StorageType* localData = dynamic_cast<StorageType*>( mLocalData.get() );

    LAMA_ASSERT_ERROR( localData, *mLocalData << ": does not fit matrix type " << typeName() )

    localData->swap( localStorage );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* COOSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<>
const char* COOSparseMatrix<float>::typeName()
{
    return "COOSparseMatrix<float>";
}

template<>
const char* COOSparseMatrix<double>::typeName()
{
    return "COOSparseMatrix<double>";
}

/* -------------------------------------------------------------------------- */
/* Template instantiation for float and double                                */
/* -------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT COOSparseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT COOSparseMatrix<double> ;

}
