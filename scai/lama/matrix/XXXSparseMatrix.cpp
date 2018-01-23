/**
 * @file XXXSparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class XXXSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 */

// hpp
#include <scai/lama/matrix/XXXSparseMatrix.hpp>

#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>

using std::shared_ptr;

namespace scai
{

using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, XXXSparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.XXXSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::shared_ptr<MatrixStorage<ValueType> > XXXSparseMatrix<ValueType>::createStorage()
{
    return shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
}

template<typename ValueType>
std::shared_ptr<MatrixStorage<ValueType> > XXXSparseMatrix<ValueType>::createStorage(
    const IndexType numRows,
    const IndexType numColumns )
{
    shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
    storage->allocate( numRows, numColumns );
    return storage;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix()

    : SparseMatrix<ValueType>( createStorage() )

{
    SCAI_LOG_INFO( logger, "XXXSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const IndexType numRows, const IndexType numColumns )

    : SparseMatrix<ValueType>( createStorage( numRows, numColumns ) )

{
    SCAI_LOG_INFO( logger, "XXXSparseMatrix( " << numRows << " x " << numColumns << " )" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

    : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ), rowDist,
                               colDist )
{
    // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const XXXSparseMatrix& other )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContextPtr( other.getContextPtr() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Matrix<ValueType>& other, bool transposeFlag ) : 

    SparseMatrix<ValueType>( createStorage() )

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
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Matrix<ValueType>& other, DistributionPtr rowDist, DistributionPtr colDist ) : 

    SparseMatrix<ValueType>( createStorage() )

{
    this->setContextPtr( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );
    // this might be done more efficiently as assign introduces intermediate copy
    SparseMatrix<ValueType>::assign( other );
    this->redistribute( rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const MatrixStorage<ValueType>& globalData ) : 

    SparseMatrix<ValueType>( createStorage() )

{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumColumns() ) );
    SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix(
    const MatrixStorage<ValueType>& localData,
    DistributionPtr rowDist,
    DistributionPtr colDist ) : 

    SparseMatrix<ValueType>( createStorage() )

{
    SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Expression_SM<ValueType>& expression ) : 

    SparseMatrix<ValueType>( createStorage() )

{
    const Matrix<ValueType>& master = expression.getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Expression_SMM<ValueType>& expression ) : 

    SparseMatrix<ValueType>( createStorage() )

{
    const Matrix<ValueType>& master = expression.getArg1().getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Expression_SM_SM<ValueType>& expression )

    : SparseMatrix<ValueType>( createStorage() )
{
    // inherit context from matA in alpha * matA + beta * matB
    const Matrix<ValueType>& master = expression.getArg1().getArg2();
    SparseMatrix<ValueType>::setContextPtr( master.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( master.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const std::string& filename )

    : SparseMatrix<ValueType>( createStorage() )

{
    this->readFromFile( filename );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::~XXXSparseMatrix()
{
    SCAI_LOG_INFO( logger, "~XXXSpareMatrix" )
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
XXXSparseMatrix<ValueType>& XXXSparseMatrix<ValueType>::operator=( const XXXSparseMatrix& matrix )
{
    SCAI_LOG_INFO( logger, "XXXSparseMatrix = XXXSparseMatrix : " << matrix )
    this->assign( matrix );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename XXXSparseMatrix<ValueType>::StorageType&
XXXSparseMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type
    const StorageType* local = dynamic_cast<const StorageType*>( this->mLocalData.get() );
    SCAI_ASSERT_ERROR( local, "XXXSparseMatrix: local storage is no more XXX: " << *this->mLocalData )
    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
typename XXXSparseMatrix<ValueType>::StorageType&
XXXSparseMatrix<ValueType>::getLocalStorage()
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type
    StorageType* local = dynamic_cast<StorageType*>( this->mLocalData.get() );
    SCAI_ASSERT_ERROR( local, "XXXSparseMatrix: local storage is no more XXX: " << *this->mLocalData )
    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const typename XXXSparseMatrix<ValueType>::StorageType&
XXXSparseMatrix<ValueType>::getHaloStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type
    const StorageType* halo = dynamic_cast<const StorageType*>( mHaloData.get() );
    SCAI_ASSERT_ERROR( halo, "XXXSparseMatrix: halo storage is no more XXX: " << *mHaloData )
    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void XXXSparseMatrix<ValueType>::swapLocalStorage( StorageType& localStorage )
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
XXXSparseMatrix<ValueType>* XXXSparseMatrix<ValueType>::newMatrix() const
{
    std::unique_ptr<XXXSparseMatrix<ValueType> > newSparseMatrix( new XXXSparseMatrix<ValueType>() );
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
XXXSparseMatrix<ValueType>* XXXSparseMatrix<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy of " << *this )
    XXXSparseMatrix<ValueType>* newSparseMatrix = new XXXSparseMatrix<ValueType>( *this );
    SCAI_LOG_INFO( logger, "copy is " << *newSparseMatrix )
    return newSparseMatrix;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* XXXSparseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
_Matrix* XXXSparseMatrix<ValueType>::create()
{
    return new XXXSparseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType XXXSparseMatrix<ValueType>::createValue()
{
    return MatrixCreateKeyType( Format::XXX, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string XXXSparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "XXXSparseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* XXXSparseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( XXXSparseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
