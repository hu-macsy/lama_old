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
#include <scai/lama/examples/stencil/OpenMPStencilKernel.hpp>

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
using namespace tasking;

using common::unique_ptr;
using common::shared_ptr;
using common::TypeTraits;
using common::binary;

namespace lama
{

template<typename ValueType>
StencilStorage<ValueType>::StencilStorage( const common::Grid& grid, const Stencil<ValueType>&  stencil ) :

    MatrixStorage<ValueType>( grid.size(), grid.size() ),
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

#define SCAI_STENCIL_MAX_POINTS 128

template<typename ValueType>
void StencilStorage<ValueType>::buildCSRData( HArray<IndexType>& csrIA, HArray<IndexType>& csrJA, _HArray& csrValues ) const
{
    IndexType n = this->getNumRows();

    int       stencilPos[SCAI_STENCIL_MAX_POINTS];

    IndexType gridDistances[SCAI_GRID_MAX_DIMENSION];

    mGrid.getDistances( gridDistances );

    mStencil.getLinearPositions( stencilPos, gridDistances );

    const IndexType* gridSizes = mGrid.sizes();

    for ( IndexType i = 0; i < mStencil.nPoints(); ++i )
    {
        SCAI_LOG_DEBUG( logger, "point = " << i << ", pos = " << stencilPos[i] << ", val = " << mStencil.values()[i] )
    }

    {
        WriteOnlyAccess<IndexType> sizes( csrIA, n );
  
        for ( IndexType i = 0; i < n; ++i )
        {
            IndexType gridPos[SCAI_GRID_MAX_DIMENSION];
            bool      valid [SCAI_STENCIL_MAX_POINTS];

            mGrid.gridPos( gridPos, i );   // grid position of point i 

            sizes[i] = mStencil.getValidPoints( valid, gridSizes, gridPos );

            SCAI_LOG_TRACE( logger, "Grid point " << i << " has " << sizes[i] << " valid neighbors" )
        }
    }

    // now build offset array

    IndexType nnz = HArrayUtils::scan1( csrIA );

    SCAI_LOG_DEBUG( logger, "Stencil -> CSR, nnz = " << nnz )

    // compute ja and values array

    HArray<ValueType> typedValues;
 
    const ValueType* stencilValues = mStencil.values();

    {
        ReadAccess<IndexType> ia( csrIA );
        WriteOnlyAccess<IndexType> ja( csrJA, nnz );
        WriteOnlyAccess<ValueType> values( typedValues, nnz );
 
        for ( IndexType i = 0; i < n; ++i )
        {
            IndexType gridPos[SCAI_GRID_MAX_DIMENSION];
            bool      valid [SCAI_STENCIL_MAX_POINTS];

            mGrid.gridPos( gridPos, i );   // grid position of point i 

            mStencil.getValidPoints( valid, gridSizes, gridPos );

            IndexType pos = ia[i];

            for ( IndexType k = 0; k < mStencil.nPoints(); ++k )
            {
                if ( valid[k] )
                {
                    ja[pos]     = i + stencilPos[k];
                    values[pos] = stencilValues[k];
                    pos++;
                }
            }
 
            SCAI_ASSERT_EQ_ERROR( pos, ia[i+1], "serious mismatch" )
        }
    }

    if ( typedValues.getValueType() == csrValues.getValueType() )
    {
        typedValues.swap( csrValues );
    }
    else
    {
        HArrayUtils::assign( csrValues, typedValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* StencilStorage<ValueType>::incGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool /* async */ ) const
{
    ReadAccess<ValueType> rX( x );
    WriteAccess<ValueType> wResult( result );

    IndexType gridDistances[SCAI_GRID_MAX_DIMENSION];

    IndexType lb[SCAI_GRID_MAX_DIMENSION];
    IndexType ub[SCAI_GRID_MAX_DIMENSION];

    mGrid.getDistances( gridDistances );

    mStencil.getWidths( lb, ub );

    int stencilLinPos[SCAI_STENCIL_MAX_POINTS];

    mStencil.getLinearPositions( stencilLinPos, gridDistances );

    stencilkernel::OpenMPStencilKernel::stencilGEMV( wResult.get(), alpha, rX.get(), 
                 mGrid.nDims(), mGrid.sizes(), lb, ub, gridDistances,
                 mStencil.nPoints(), mStencil.positions(), mStencil.values(),
                 stencilLinPos );

    return NULL;  
}

template<typename ValueType>
void StencilStorage<ValueType>::matrixTimesVector(

    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const

{
    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x 
                         << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.Stencil.matrixTimesVector" )

    if ( alpha == common::constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, common::binary::MULT, y, this->getContextPtr() );
        return;
    }

    ContextPtr loc = this->getContextPtr();

    // GEMV only implemented as y += A * x, so split

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumRows );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, common::binary::MULT, y, this->getContextPtr() );
    }

    bool async = false;

    SyncToken* token = incGEMV( result, alpha, x, async );

    SCAI_ASSERT( token == NULL, "syncrhonous execution cannot have token" )
}

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, StencilStorage<ValueType>::logger, "MatrixStorage.StencilStorage" )

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StencilStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
