/**
 * @file JDSSparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of methods and constructors for template class JDSSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 */

// hpp
#include <scai/lama/matrix/JDSSparseMatrix.hpp>

#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using common::shared_ptr;
using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSSparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.JDSSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage(
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
    SCAI_LOG_INFO( logger, "JDSSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    SCAI_LOG_INFO( logger, "JDSSparseMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ), rowDist,
                               colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const JDSSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContextPtr( other.getContextPtr() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Matrix& other, bool transposeFlag )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setContextPtr( other.getContextPtr() );
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
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setContextPtr( other.getContextPtr() );
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

    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );

    Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Expression_SMM& expression )

    : SparseMatrix<ValueType>( createStorage() )

{
    const Matrix& master = expression.getArg1().getArg2();

    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
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

    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
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
    SCAI_LOG_INFO( logger, "~JDSSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
JDSSparseMatrix<ValueType>& JDSSparseMatrix<ValueType>::operator=( const JDSSparseMatrix& matrix )
{
    SCAI_LOG_INFO( logger, "JDSSparseMatrix = JDSSparseMatrix : " << matrix )
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

    SCAI_ASSERT_ERROR( local, "JDSSparseMatrix: local storage is no more JDS: " << *this->mLocalData )

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

    SCAI_ASSERT_ERROR( local, "JDSSparseMatrix: local storage is no more JDS: " << *this->mLocalData )

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

    SCAI_ASSERT_ERROR( halo, "JDSSparseMatrix: halo storage is no more JDS: " << *mHaloData )

    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
{
    // make sure that local storage fits into this sparse matrix

    SCAI_ASSERT_EQUAL_ERROR( localStorage.getNumRows(), mLocalData->getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( localStorage.getNumColumns(), mLocalData->getNumColumns() )

    // make sure that local matrix storage has the correct format / value type

    StorageType* localData = dynamic_cast<StorageType*>( mLocalData.get() );

    SCAI_ASSERT_ERROR( localData, *mLocalData << ": does not fit matrix type " << typeName() )

    localData->swap( localStorage );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::newMatrix() const
{
    common::unique_ptr<JDSSparseMatrix<ValueType> > newSparseMatrix( new JDSSparseMatrix<ValueType>() );

    // inherit the context, communication kind of this matrix for the new matrix

    newSparseMatrix->setContextPtr( this->getContextPtr() );

    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );

    SCAI_LOG_INFO( logger,
                   *this << ": create -> " << *newSparseMatrix << " @ " << *(newSparseMatrix->getContextPtr())
                   << ", kind = " << newSparseMatrix->getCommunicationKind() );

    return newSparseMatrix.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy of " << *this )

    JDSSparseMatrix<ValueType>* newSparseMatrix = new JDSSparseMatrix<ValueType>( *this );

    SCAI_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* JDSSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Matrix* JDSSparseMatrix<ValueType>::create()
{
    return new JDSSparseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType JDSSparseMatrix<ValueType>::createValue()
{
    return MatrixCreateKeyType( Format::JDS, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string JDSSparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("JDSSparseMatrix<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* JDSSparseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( JDSSparseMatrix, SCAI_ARITHMETIC_HOST )

} /* end namespace lama */

} /* end namespace scai */
