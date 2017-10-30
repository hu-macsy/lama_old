/**
 * @file COOSparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class COOSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 */

// hpp
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using common::shared_ptr;
using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOSparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.COOSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > COOSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
common::shared_ptr<MatrixStorage<ValueType> > COOSparseMatrix<ValueType>::createStorage(
    const IndexType numRows,
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    SCAI_LOG_INFO( logger, "COOSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    SCAI_LOG_INFO( logger, "COOSparseMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ), rowDist,
                               colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const COOSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContextPtr( other.getContextPtr() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const _Matrix& other, bool transposeFlag )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setContextPtr( other.getContextPtr() );
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
COOSparseMatrix<ValueType>::COOSparseMatrix( const _Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )

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
COOSparseMatrix<ValueType>::COOSparseMatrix( const _MatrixStorage& globalData )

    : SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumColumns() ) );
    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix(
    const _MatrixStorage& localData,
    DistributionPtr rowDist,
    DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const Expression_SM& expression )

    : SparseMatrix<ValueType>( createStorage() )

{
    const _Matrix& master = expression.getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    _Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const Expression_SMM& expression )

    : SparseMatrix<ValueType>( createStorage() )

{
    const _Matrix& master = expression.getArg1().getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    _Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const Expression_SM_SM& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB
    const _Matrix& master = expression.getArg1().getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    _Matrix::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::~COOSparseMatrix()
{
    SCAI_LOG_INFO( logger, "~COOSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
COOSparseMatrix<ValueType>& COOSparseMatrix<ValueType>::operator=( const COOSparseMatrix& matrix )
{
    SCAI_LOG_INFO( logger, "COOSparseMatrix = COOSparseMatrix : " << matrix )
    this->assign( matrix );
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
    SCAI_ASSERT_ERROR( local, "COOSparseMatrix: local storage is no more COO: " << *this->mLocalData )
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
    SCAI_ASSERT_ERROR( local, "COOSparseMatrix: local storage is no more COO: " << *this->mLocalData )
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
    SCAI_ASSERT_ERROR( halo, "COOSparseMatrix: halo storage is no more COO: " << *mHaloData )
    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
COOSparseMatrix<ValueType>* COOSparseMatrix<ValueType>::newMatrix() const
{
    common::unique_ptr<COOSparseMatrix<ValueType> > newSparseMatrix( new COOSparseMatrix<ValueType>() );
    // inherit the context, communication kind of this matrix for the new matrix
    newSparseMatrix->setContextPtr( this->getContextPtr() );
    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );
    SCAI_LOG_INFO( logger,
                   *this << ": create -> " << *newSparseMatrix << " @ " << * ( newSparseMatrix->getContextPtr() )
                   << ", kind = " << newSparseMatrix->getCommunicationKind() );
    return newSparseMatrix.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>* COOSparseMatrix<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy of " << *this )
    COOSparseMatrix<ValueType>* newSparseMatrix = new COOSparseMatrix<ValueType>( *this );
    SCAI_LOG_INFO( logger, "copy is " << *newSparseMatrix )
    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* COOSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
_Matrix* COOSparseMatrix<ValueType>::create()
{
    return new COOSparseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType COOSparseMatrix<ValueType>::createValue()
{
    return MatrixCreateKeyType( Format::COO, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string COOSparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "COOSparseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* COOSparseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( COOSparseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
