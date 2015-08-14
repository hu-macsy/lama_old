/**
 * @file DIASparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 */

// hpp
#include <scai/lama/matrix/DIASparseMatrix.hpp>

// boost
#include <boost/preprocessor.hpp>

using common::shared_ptr;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIASparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.DIASparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage(
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
    SCAI_LOG_INFO( logger, "DIASparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    SCAI_LOG_INFO( logger, "DIASparseMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ), rowDist,
                               colDist )
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
    this->setContext( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    if( transposeFlag )
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
DIASparseMatrix<ValueType>::DIASparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )

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
DIASparseMatrix<ValueType>::DIASparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumColumns() ) );

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
    const Matrix& master = expression.getArg2();

    SparseMatrix<ValueType>::setContext( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Expression_SMM& expression )

    : SparseMatrix<ValueType>( createStorage() )

{
    const Matrix& master = expression.getArg1().getArg2();

    SparseMatrix<ValueType>::setContext( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Expression_SM_SM& expression )

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
DIASparseMatrix<ValueType>::DIASparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::~DIASparseMatrix()
{
    SCAI_LOG_INFO( logger, "~DIASpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
DIASparseMatrix<ValueType>& DIASparseMatrix<ValueType>::operator=( const DIASparseMatrix& matrix )
{
    SCAI_LOG_INFO( logger, "DIASparseMatrix = DIASparseMatrix : " << matrix )
    this->assign( matrix );
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
DIASparseMatrix<ValueType>* DIASparseMatrix<ValueType>::clone() const
{
    DIASparseMatrix* newSparseMatrix = new DIASparseMatrix<ValueType>();

    // inherit the context, communication kind of this matrix for the new matrix

    newSparseMatrix->setContext( this->getContextPtr() );
    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );

    SCAI_LOG_INFO( logger, "create is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>* DIASparseMatrix<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy of " << *this )

    DIASparseMatrix<ValueType>* newSparseMatrix = new DIASparseMatrix<ValueType>( *this );

    SCAI_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* DIASparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Matrix* DIASparseMatrix<ValueType>::create()
{
    return new DIASparseMatrix<ValueType>();
}

template<typename ValueType>
std::pair<MatrixStorageFormat, common::ScalarType> DIASparseMatrix<ValueType>::createValue()
{
    common::ScalarType skind = common::getScalarType<ValueType>();
    return std::pair<MatrixStorageFormat, common::ScalarType> ( Format::DIA, skind );
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

#define tmp_xstr(s) tmp_str(s)
#define tmp_str(s) #s

#define LAMA_DIA_SPARSE_MATRIX_INSTANTIATE(z, I, _)                             \
                                                                                \
    template<>                                                                  \
    const char* DIASparseMatrix<ARITHMETIC_TYPE##I>::typeName()                 \
    {                                                                           \
        return "DIASparseMatrix<" tmp_xstr(ARITHMETIC_TYPE##I) ">";             \
    }                                                                           \
                                                                                \
    template class COMMON_DLL_IMPORTEXPORT DIASparseMatrix<ARITHMETIC_TYPE##I> ;

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_DIA_SPARSE_MATRIX_INSTANTIATE, _ )

#undef LAMA_DIA_SPARSE_MATRIX_INSTANTIATE
#undef tmp_xstr
#undef tmp_str

}
