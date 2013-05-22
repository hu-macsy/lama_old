/**
 * @file DIASparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class DIASparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 * $Id$
 */

// hpp
#include <lama/matrix/DIASparseMatrix.hpp>

using boost::shared_ptr;

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, DIASparseMatrix<T>::logger, "Matrix.SparseMatrix.DIASparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage(
    const IndexType numRows, 
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    LAMA_LOG_INFO( logger, "DIASpareMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    LAMA_LOG_INFO( logger, "DIASpareMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ),
                                   rowDist, colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const DIASparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContext( other.getContextPtr() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Matrix& other, bool transposeFlag )

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
DIASparseMatrix<ValueType>::DIASparseMatrix(
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
DIASparseMatrix<ValueType>::DIASparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumRows() ) );

    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix(
    const _MatrixStorage& localData, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Expression_SM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Expression_SMM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Expression_SM_SM& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB

    SparseMatrix<ValueType>::setContext( expression.getArg1().getArg2().getContextPtr() );
    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::~DIASparseMatrix()
{
    LAMA_LOG_INFO( logger, "~DIASpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
DIASparseMatrix<ValueType>& DIASparseMatrix<ValueType>::operator=( const DIASparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, "DIASparseMatrix = DIASparseMatrix : " << matrix )
    assign( matrix );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename DIASparseMatrix<ValueType>::StorageType&
DIASparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "DIASparseMatrix: local storage is no more DIA: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename DIASparseMatrix<ValueType>::StorageType&
DIASparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "DIASparseMatrix: local storage is no more DIA: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename DIASparseMatrix<ValueType>::StorageType&
DIASparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );

    LAMA_ASSERT_ERROR( halo, "DIASparseMatrix: halo storage is no more DIA: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIASparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
DIASparseMatrix<ValueType>* DIASparseMatrix<ValueType>::create() const
{
    DIASparseMatrix<ValueType>* newSparseMatrix = new DIASparseMatrix<ValueType>();

    // inherit the context of this matrix for the new matrix

    newSparseMatrix->setContext( this->getContextPtr() );

    LAMA_LOG_INFO( logger, "create is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>* DIASparseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy of " << *this )

    DIASparseMatrix<ValueType>* newSparseMatrix = new DIASparseMatrix<ValueType>( *this );

    LAMA_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* DIASparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<>
const char* DIASparseMatrix<float>::typeName()
{
    return "DIASparseMatrix<float>";
}

template<>
const char* DIASparseMatrix<double>::typeName()
{
    return "DIASparseMatrix<double>";
}

/* -------------------------------------------------------------------------- */
/* Template instantiation for float and double                                */
/* -------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT DIASparseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT DIASparseMatrix<double> ;

}
