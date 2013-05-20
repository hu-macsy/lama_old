/**
 * @file CSRSparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class CSRSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 * $Id$
 */

// hpp
#include <lama/matrix/CSRSparseMatrix.hpp>

using boost::shared_ptr;

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, CSRSparseMatrix<T>::logger, "Matrix.SparseMatrix.CSRSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > CSRSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > CSRSparseMatrix<ValueType>::createStorage(
    const IndexType numRows, 
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    LAMA_LOG_INFO( logger, "CSRSpareMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    LAMA_LOG_INFO( logger, "CSRSpareMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ),
                                   rowDist, colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const CSRSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContext( other.getContextPtr() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const Matrix& other, bool transposeFlag )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );

    if ( transposeFlag )
    {
        SparseMatrix<ValueType>::assignTranspose( other );
    }
    else
    {
        SparseMatrix<ValueType>::assign( other );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix(
    const Matrix& other, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );

    // this might be done more efficiently as assign introduces intermediate copy

    SparseMatrix<ValueType>::assign( other );
    this->redistribute( rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumRows() ) );

    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix(
    const _MatrixStorage& localData, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const Expression<Matrix, Matrix, Times>& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( const Expression<Scalar, Matrix, Times>& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( 
    const Expression<Scalar, Expression<Matrix, Matrix, Times>, Times>& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix( 
    const Expression<Expression<Scalar, Matrix, Times>,
                     Expression<Scalar, Matrix, Times>,
                     Plus> expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB

    SparseMatrix<ValueType>::setContext( expression.getArg1().getArg2().getContextPtr() );
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::CSRSparseMatrix(const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>::~CSRSparseMatrix()
{
    LAMA_LOG_INFO( logger, "~CSRSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
CSRSparseMatrix<ValueType>& CSRSparseMatrix<ValueType>::operator=( const CSRSparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, "CSRSparseMatrix = CSRSparseMatrix : " << matrix )
    assign( matrix );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename CSRSparseMatrix<ValueType>::StorageType&
CSRSparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "CSRSparseMatrix: local storage is no more CSR: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename CSRSparseMatrix<ValueType>::StorageType&
CSRSparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "CSRSparseMatrix: local storage is no more CSR: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename CSRSparseMatrix<ValueType>::StorageType&
CSRSparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );

    LAMA_ASSERT_ERROR( halo, "CSRSparseMatrix: halo storage is no more CSR: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
CSRSparseMatrix<ValueType>* CSRSparseMatrix<ValueType>::create() const
{
    CSRSparseMatrix<ValueType>* newSparseMatrix = new CSRSparseMatrix<ValueType>();

    // inherit the context of this matrix for the new matrix

    newSparseMatrix->setContext( this->getContextPtr() );

    LAMA_LOG_INFO( logger, "create is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
CSRSparseMatrix<ValueType>* CSRSparseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy of " << *this )

    CSRSparseMatrix<ValueType>* newSparseMatrix = new CSRSparseMatrix<ValueType>( *this );

    LAMA_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* CSRSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<>
const char* CSRSparseMatrix<float>::typeName()
{
    return "CSRSparseMatrix<float>";
}

template<>
const char* CSRSparseMatrix<double>::typeName()
{
    return "CSRSparseMatrix<double>";
}

/* -------------------------------------------------------------------------- */
/* Template instantiation for float and double                                */
/* -------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT CSRSparseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT CSRSparseMatrix<double> ;

}
