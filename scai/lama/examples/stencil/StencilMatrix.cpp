/**
 * @file StencilMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class StencilMatrix.
 * @author Thomas Brandes
 * @date 13.04.2017
 */

// hpp
#include "StencilMatrix.hpp"

#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/shared_ptr.hpp>

namespace scai
{

using common::shared_ptr;
using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, StencilMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.StencilMatrix" )

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( const common::Grid& grid, const Stencil<ValueType>& stencil )

    : SparseMatrix<ValueType>( common::shared_ptr<MatrixStorage<ValueType> >( new StencilStorage<ValueType>( grid, stencil ) ) )

{
}

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( dmemo::DistributionPtr dist, const Stencil<ValueType>& stencil )

    : SparseMatrix<ValueType>( common::shared_ptr<MatrixStorage<ValueType> >( new StencilStorage<ValueType>( common::Grid( 0 ), stencil ) ) )

{
    SCAI_ASSERT_ERROR( dist, "NULL dist" )
}

template<typename ValueType>
StencilMatrix<ValueType>::~StencilMatrix()
{
    SCAI_LOG_INFO( logger, "~StencilMatrix" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const StencilStorage<ValueType>& StencilMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type
    const StencilStorage<ValueType>* local = dynamic_cast<const StencilStorage<ValueType>*>( this->mLocalData.get() );
    SCAI_ASSERT_ERROR( local, "StencilMatrix: local storage is no more Stencil: " << *this->mLocalData )
    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const CSRStorage<ValueType>& StencilMatrix<ValueType>::getHaloStorage() const
{
    const CSRStorage<ValueType>* halo = dynamic_cast<const CSRStorage<ValueType>*>( this->mHaloData.get() );
    SCAI_ASSERT_ERROR( halo, "StencilMatrix: local storage is no more CSR: " << *this->mHaloData )
    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>* StencilMatrix<ValueType>::newMatrix() const
{
    return NULL;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>* StencilMatrix<ValueType>::copy() const
{
    return NULL;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* StencilMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::string StencilMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "StencilMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* StencilMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StencilMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
