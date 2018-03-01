/**
 * @file DIASparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class DIASparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 */

// hpp
#include <scai/lama/matrix/DIASparseMatrix.hpp>

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

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIASparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.DIASparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::shared_ptr<DIAStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage( hmemo::ContextPtr ctx )
{
    return shared_ptr<DIAStorage<ValueType> >( new StorageType( ctx ) );
}

template<typename ValueType>
std::shared_ptr<DIAStorage<ValueType> > DIASparseMatrix<ValueType>::createStorage( DIAStorage<ValueType>&& storage )
{
    return shared_ptr<DIAStorage<ValueType> >( new DIAStorage<ValueType>( std::move( storage ) ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( hmemo::ContextPtr ctx ) : 

    SparseMatrix<ValueType>( createStorage( ctx ) )

{
    SCAI_LOG_INFO( logger, "DIASparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const DIASparseMatrix& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( DIASparseMatrix&& other ) noexcept :

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>&  DIASparseMatrix<ValueType>::operator=( DIASparseMatrix&& other ) 
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const Matrix<ValueType>& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setContextPtr( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( DIAStorage<ValueType> globalStorage ) : 

    SparseMatrix<ValueType>( createStorage( std::move( globalStorage ) ) )

{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( DistributionPtr rowDist, DIAStorage<ValueType> localStorage ) :

    SparseMatrix<ValueType>( rowDist, createStorage( std::move( localStorage ) ) )
{
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
    SCAI_ASSERT_ERROR( local, "DIASparseMatrix: local storage is no more DIA: " << *this->mLocalData )
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
    SCAI_ASSERT_ERROR( local, "DIASparseMatrix: local storage is no more DIA: " << *this->mLocalData )
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
    SCAI_ASSERT_ERROR( halo, "DIASparseMatrix: halo storage is no more DIA: " << *mHaloData )
    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DIASparseMatrix<ValueType>* DIASparseMatrix<ValueType>::newMatrix() const
{
    std::unique_ptr<DIASparseMatrix<ValueType> > newSparseMatrix( new DIASparseMatrix<ValueType>() );
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
_Matrix* DIASparseMatrix<ValueType>::create()
{
    return new DIASparseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType DIASparseMatrix<ValueType>::createValue()
{
    return MatrixCreateKeyType( Format::DIA, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string DIASparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "DIASparseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* DIASparseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DIASparseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
