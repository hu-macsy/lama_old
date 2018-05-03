/**
 * @file JDSSparseMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class JDSSparseMatrix.
 * @author Thomas Brandes
 * @date 04.08.2012
 */

// hpp
#include <scai/lama/matrix/JDSSparseMatrix.hpp>

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

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSSparseMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.JDSSparseMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::shared_ptr<JDSStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage( hmemo::ContextPtr ctx )
{
    return shared_ptr<JDSStorage<ValueType> >( new StorageType( ctx ) );
}

template<typename ValueType>
std::shared_ptr<JDSStorage<ValueType> > JDSSparseMatrix<ValueType>::createStorage( JDSStorage<ValueType>&& storage )
{
    return shared_ptr<JDSStorage<ValueType> >( new JDSStorage<ValueType>( std::move( storage ) ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( hmemo::ContextPtr ctx ) : 

    SparseMatrix<ValueType>( createStorage( ctx ) )

{
    SCAI_LOG_INFO( logger, "JDSSparseMatrix()" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const JDSSparseMatrix& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setCommunicationKind( other.getCommunicationKind() );
    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( JDSSparseMatrix&& other ) noexcept :

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>&  JDSSparseMatrix<ValueType>::operator=( JDSSparseMatrix&& other ) 
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( const Matrix<ValueType>& other ) : 

    SparseMatrix<ValueType>( createStorage( other.getContextPtr() ) )

{
    this->setContextPtr( other.getContextPtr() );
    this->setCommunicationKind( other.getCommunicationKind() );

    SparseMatrix<ValueType>::assign( other );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( JDSStorage<ValueType> globalStorage ) : 

    SparseMatrix<ValueType>( createStorage( std::move( globalStorage ) ) )

{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
JDSSparseMatrix<ValueType>::JDSSparseMatrix( DistributionPtr rowDist, JDSStorage<ValueType> localStorage ) :

    SparseMatrix<ValueType>( rowDist, createStorage( std::move( localStorage ) ) )
{
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
JDSSparseMatrix<ValueType>* JDSSparseMatrix<ValueType>::newMatrix() const
{
    std::unique_ptr<JDSSparseMatrix<ValueType> > newSparseMatrix( new JDSSparseMatrix<ValueType>() );
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
_Matrix* JDSSparseMatrix<ValueType>::create()
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
    s << std::string( "JDSSparseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
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

SCAI_COMMON_INST_CLASS( JDSSparseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
