/**
 * @file GridVector.cpp
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
 * @brief Implementation of methods for grid vector
 * @author Thomas Brandes
 * @date 12.05.2017
 */

#include <scai/lama/GridVector.hpp>

// other SCAI libraries

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

template<typename ValueType>
void GridVector<ValueType>::reduce( const GridVector<ValueType>& other, IndexType dim, const common::binary::BinaryOp redOp )
{
    SCAI_ASSERT_VALID_INDEX_ERROR( dim, other.nDims(), "illegeal reduction dim on this grid " << other.globalGrid() )

    COMMON_THROWEXCEPTION( "reduction on grids not supported yet, reduction op = " << redOp )
}

template<typename ValueType>
void GridVector<ValueType>::gemm( const ValueType alpha, const GridVector<ValueType>& v1, const GridVector<ValueType>& v2 )
{
    const common::Grid& resGrid = this->globalGrid();
    const common::Grid& grid1 = v1.globalGrid();
    const common::Grid& grid2 = v2.globalGrid();

    SCAI_ASSERT_EQ_ERROR( 2, resGrid.nDims(), "gemm only on two-dimensional grids." )
    SCAI_ASSERT_EQ_ERROR( 2, grid1.nDims(), "gemm only on two-dimensional grids." )
    SCAI_ASSERT_EQ_ERROR( 2, grid2.nDims(), "gemm only on two-dimensional grids." )

    SCAI_ASSERT_EQ_ERROR( resGrid.size( 0 ), grid1.size( 0 ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( resGrid.size( 1 ), grid2.size( 1 ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( grid1.size( 1 ), grid2.size( 0 ), "size mismatch" )

    // ToDo: not yet for distributed grids

    const IndexType n = resGrid.size( 0 );
    const IndexType m = resGrid.size( 1 );
    const IndexType k = grid1.size( 0 );

    int lda = grid1.size( 1 );
    int ldb = grid2.size( 1 );
    int ldc = resGrid.size( 1 );

    if ( lda != 0 && n != 0 && m != 0 )
    {
        static utilskernel::LAMAKernel<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;

        ContextPtr loc = Context::getHostPtr();
        gemm.getSupportedContext( loc );

        GridWriteAccess<ValueType> wRes( *this, loc );
        GridReadAccess<ValueType> rA( v1, loc );
        GridReadAccess<ValueType> rB( v2, loc );
  
        SCAI_CONTEXT_ACCESS( loc )

        ValueType beta = 1;

        gemm[loc]( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, rA.get(), lda, rB.get(), ldb, beta,
                   wRes.get(), ldc );
    }
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( GridVector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
