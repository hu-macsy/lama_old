/**
 * @file XXXSparseMatrix.cpp
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
std::shared_ptr<XXXStorage<ValueType> > XXXSparseMatrix<ValueType>::createStorage( hmemo::ContextPtr ctx )
{
    return shared_ptr<XXXStorage<ValueType> >( new StorageType( ctx ) );
}

template<typename ValueType>
std::shared_ptr<XXXStorage<ValueType> > XXXSparseMatrix<ValueType>::createStorage( XXXStorage<ValueType>&& storage )
{
    return shared_ptr<XXXStorage<ValueType> >( new XXXStorage<ValueType>( std::move( storage ) ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( hmemo::ContextPtr ctx ) : 

    SparseMatrix<ValueType>( createStorage( ctx ) )

{
    SCAI_LOG_INFO( logger, "XXXSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const XXXSparseMatrix& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( XXXSparseMatrix&& other ) noexcept :

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>&  XXXSparseMatrix<ValueType>::operator=( XXXSparseMatrix&& other ) 
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( const Matrix<ValueType>& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setContextPtr( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( XXXStorage<ValueType> globalStorage ) : 

    SparseMatrix<ValueType>( createStorage( std::move( globalStorage ) ) )

{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
XXXSparseMatrix<ValueType>::XXXSparseMatrix( DistributionPtr rowDist, XXXStorage<ValueType> localStorage ) :

    SparseMatrix<ValueType>( rowDist, createStorage( std::move( localStorage ) ) )
{
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
XXXSparseMatrix<ValueType>* XXXSparseMatrix<ValueType>::newMatrix() const
{
    std::unique_ptr<XXXSparseMatrix<ValueType> > newSparseMatrix( new XXXSparseMatrix<ValueType>() );
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
