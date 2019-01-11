/**
 * @file COOSparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

#include <memory>

using std::shared_ptr;

namespace scai
{

using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOSparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.COOSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::shared_ptr<COOStorage<ValueType> > COOSparseMatrix<ValueType>::createStorage( hmemo::ContextPtr ctx )
{
    return shared_ptr<COOStorage<ValueType> >( new StorageType( ctx ) );
}

template<typename ValueType>
std::shared_ptr<COOStorage<ValueType> > COOSparseMatrix<ValueType>::createStorage( COOStorage<ValueType>&& storage )
{
    return shared_ptr<COOStorage<ValueType> >( new COOStorage<ValueType>( std::move( storage ) ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( hmemo::ContextPtr ctx ) : 

    SparseMatrix<ValueType>( createStorage( ctx ) )

{
    SCAI_LOG_INFO( logger, "COOSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const COOSparseMatrix& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( COOSparseMatrix&& other ) noexcept :

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>&  COOSparseMatrix<ValueType>::operator=( COOSparseMatrix&& other ) 
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( const Matrix<ValueType>& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setContextPtr( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( COOStorage<ValueType> globalStorage ) : 

    SparseMatrix<ValueType>( createStorage( std::move( globalStorage ) ) )

{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix( DistributionPtr rowDist, COOStorage<ValueType> localStorage ) :

    SparseMatrix<ValueType>( rowDist, createStorage( std::move( localStorage ) ) )
{
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
COOSparseMatrix<ValueType>* COOSparseMatrix<ValueType>::newMatrix() const
{
    std::unique_ptr<COOSparseMatrix<ValueType> > newSparseMatrix( new COOSparseMatrix<ValueType>() );
    // inherit the context, communication kind of this matrix for the new matrix
    newSparseMatrix->setContextPtr( this->getContextPtr() );
    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );
    newSparseMatrix->allocate( getRowDistributionPtr(), getColDistributionPtr() );
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
