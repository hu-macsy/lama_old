/**
 * @file StencilStorage.cpp
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
 * @brief Implementation and instantiation for template class StencilStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 */

// hpp
#include "StencilStorage.hpp"

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>

#include <scai/lama/storage/StorageMethods.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/dmemo/Redistributor.hpp>


// internal scai libraries
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;
using namespace dmemo;
using namespace utilskernel;

using common::unique_ptr;
using common::shared_ptr;
using common::TypeTraits;
using common::binary;

namespace lama
{

template<typename ValueType>
StencilStorage<ValueType>::StencilStorage( const common::Grid& grid, const Stencil<ValueType>&  stencil ) :

    MatrixStorage<ValueType>(),
    mGrid( grid ),
    mStencil( stencil )
{
}

template<typename ValueType>
StencilStorage<ValueType>::~StencilStorage()
{
}

template<typename ValueType>
MatrixStorageCreateKeyType StencilStorage<ValueType>::getCreateValue() const
{
    return MatrixStorageCreateKeyType( Format::CSR, common::getScalarType<ValueType>() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::conj()
{
    COMMON_THROWEXCEPTION( "conj unsupported" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType StencilStorage<ValueType>::l1Norm() const
{
    COMMON_THROWEXCEPTION( "l1Norm unsupported" )
    return 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType StencilStorage<ValueType>::l2Norm() const
{
    COMMON_THROWEXCEPTION( "l2Norm unsupported" )
    return 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType StencilStorage<ValueType>::maxNorm() const
{
    COMMON_THROWEXCEPTION( "maxNorm unsupported" )
    return 0;
}

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, StencilStorage<ValueType>::logger, "MatrixStorage.StencilStorage" )

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StencilStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
