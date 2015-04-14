/**
 * @file JDSSparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class JDSSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 * @since 1.0.0
 */

// hpp
#include <lama/matrix/JDSSparseMatrix.hpp>

// boost
#include <boost/preprocessor.hpp>

using boost::shared_ptr;

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSSparseMatrix<ValueType>::logger, "Matrix.SparseMatrix.JDSSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
boost::shared_ptr<MatrixStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage(
    const IndexType numRows, 
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    LAMA_LOG_INFO( logger, "JDSSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    LAMA_LOG_INFO( logger, "JDSSparseMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ),
                                   rowDist, colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const JDSSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContext( other.getContextPtr() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Matrix& other, bool transposeFlag )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setContext( other.getContextPtr() );
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
JDSSparseMatrix<ValueType>::JDSSparseMatrix(
    const Matrix& other, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setContext( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    // this might be done more efficiently as assign introduces intermediate copy

    SparseMatrix<ValueType>::assign( other );
    this->redistribute( rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumColumns() ) );

    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix(
    const _MatrixStorage& localData, 
    DistributionPtr rowDist, 
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Expression_SM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    const Matrix& master = expression.getArg2(); 

    SparseMatrix<ValueType>::setContext( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Expression_SMM& expression ) 

    : SparseMatrix<ValueType>( createStorage() )

{
    const Matrix& master = expression.getArg1().getArg2(); 

    SparseMatrix<ValueType>::setContext( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Expression_SM_SM& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB

    const Matrix& master = expression.getArg1().getArg2();

    SparseMatrix<ValueType>::setContext( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::~JDSSparseMatrix()
{
    LAMA_LOG_INFO( logger, "~JDSSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
JDSSparseMatrix<ValueType>& JDSSparseMatrix<ValueType>::operator=( const JDSSparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, "JDSSparseMatrix = JDSSparseMatrix : " << matrix )
    this->assign( matrix );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename JDSSparseMatrix<ValueType>::StorageType&
JDSSparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "JDSSparseMatrix: local storage is no more JDS: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename JDSSparseMatrix<ValueType>::StorageType&
JDSSparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );

    LAMA_ASSERT_ERROR( local, "JDSSparseMatrix: local storage is no more JDS: " << *this->mLocalData )

    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename JDSSparseMatrix<ValueType>::StorageType&
JDSSparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type

    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );

    LAMA_ASSERT_ERROR( halo, "JDSSparseMatrix: halo storage is no more JDS: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::createMatrix()
{
    JDSSparseMatrix<ValueType>* newSparseMatrix = new JDSSparseMatrix<ValueType>();
    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::create() const
{
    JDSSparseMatrix* newSparseMatrix = createMatrix();

    // inherit the context, communication kind of this matrix for the new matrix

    newSparseMatrix->setContext( this->getContextPtr() );
    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );

    LAMA_LOG_INFO( logger, "create is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy of " << *this )

    JDSSparseMatrix<ValueType>* newSparseMatrix = new JDSSparseMatrix<ValueType>( *this );

    LAMA_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* JDSSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

template<typename ValueType>
bool JDSSparseMatrix<ValueType>::registerCreator()
{
    MatrixStorageFormat storageFormat = Format::JDS;

    // conversion needed even if createMatrix has only covariant return type 

    Matrix::CreateFn create = ( Matrix::CreateFn ) ( &JDSSparseMatrix<ValueType>::createMatrix );

    Matrix::addCreator( storageFormat, Scalar::getType<ValueType>(), create );

    return true;
}

template<typename ValueType>
bool JDSSparseMatrix<ValueType>::initialized = registerCreator();

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

#define LAMA_JDS_SPARSE_MATRIX_INSTANTIATE(z, I, _)                        \
                                                                           \
template<>                                                                 \
const char* JDSSparseMatrix<ARITHMETIC_TYPE##I>::typeName()                \
{                                                                          \
    return "JDSSparseMatrix<ARITHMETIC_TYPE##I>";                          \
}                                                                          \
                                                                           \
template class LAMA_DLL_IMPORTEXPORT JDSSparseMatrix<ARITHMETIC_TYPE##I> ;  

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_JDS_SPARSE_MATRIX_INSTANTIATE, _ )


#undef LAMA_JDS_SPARSE_MATRIX_INSTANTIATE

}
