/**
 * @file ELLSparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation of methods and constructors for template class ELLSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 * @since 1.0.0
 */

// hpp
#include <lama/matrix/ELLSparseMatrix.hpp>

using boost::shared_ptr;

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, ELLSparseMatrix<T>::logger, "Matrix.SparseMatrix.ELLSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > ELLSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > ELLSparseMatrix<ValueType>::createStorage(
    const IndexType numRows, 
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    LAMA_LOG_INFO( logger, "ELLSpareMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    LAMA_LOG_INFO( logger, "ELLSpareMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ),
                                   rowDist, colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const ELLSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContext( other.getContextPtr() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const Matrix& other, bool transposeFlag )

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
ELLSparseMatrix<ValueType>::ELLSparseMatrix(
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
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumRows() ) );

    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix(
    const _MatrixStorage& localData, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const Expression_SM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const Expression_SMM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const Expression_SM_SM& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB

    SparseMatrix<ValueType>::setContext( expression.getArg1().getArg2().getContextPtr() );
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>::~ELLSparseMatrix()
{
    LAMA_LOG_INFO( logger, "~ELLSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ELLSparseMatrix<ValueType>& ELLSparseMatrix<ValueType>::operator=( const ELLSparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, "ELLSparseMatrix = ELLSparseMatrix : " << matrix )
    this->assign( matrix );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename ELLSparseMatrix<ValueType>::StorageType&
ELLSparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "ELLSparseMatrix: local storage is no more ELL: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename ELLSparseMatrix<ValueType>::StorageType&
ELLSparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "ELLSparseMatrix: local storage is no more ELL: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename ELLSparseMatrix<ValueType>::StorageType&
ELLSparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );

    LAMA_ASSERT_ERROR( halo, "ELLSparseMatrix: halo storage is no more ELL: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
ELLSparseMatrix<ValueType>* ELLSparseMatrix<ValueType>::create() const
{
    ELLSparseMatrix<ValueType>* newSparseMatrix = new ELLSparseMatrix<ValueType>();

    // inherit the context of this matrix for the new matrix

    newSparseMatrix->setContext( this->getContextPtr() );

    LAMA_LOG_INFO( logger, "create is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ELLSparseMatrix<ValueType>* ELLSparseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy of " << *this )

    ELLSparseMatrix<ValueType>* newSparseMatrix = new ELLSparseMatrix<ValueType>( *this );

    LAMA_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* ELLSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<>
const char* ELLSparseMatrix<float>::typeName()
{
    return "ELLSparseMatrix<float>";
}

template<>
const char* ELLSparseMatrix<double>::typeName()
{
    return "ELLSparseMatrix<double>";
}

/* -------------------------------------------------------------------------- */
/* Template instantiation for float and double                                */
/* -------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT ELLSparseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT ELLSparseMatrix<double> ;

}
