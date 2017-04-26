/**
 * @file OpenMPStencilKernel.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief OpenMP Implementations on host for stencil kernels.
 * @author Thomas Brandes
 * @date Apr 28, 2017
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/common/OpenMP.hpp>

#include <scai/lama/examples/stencil/OpenMPStencilKernel.hpp>
#include <scai/lama/examples/stencil/StencilKernelTrait.hpp>

namespace scai
{
namespace stencilkernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPStencilKernel::logger, "OpenMP.StencilKernel" )

/* --------------------------------------------------------------------------- */

/** This routine checks whether a point pos is an inner point 
 *
 *  @param[in] pos is the position, in range 0 .. size-1
 *  @param[in] border specfies a border, if < 0 from left, if > 0 at right
 *  @param[in] size is the size of the range, only needed for border > 0
 */
static inline bool isInner( const IndexType pos, const int border, const IndexType size )
{
    if ( border == 0 )
    {
        return true;
    }
    if ( border < 0 )
    {
        return pos >= static_cast<IndexType>( -border );
    }

    return pos + static_cast<IndexType>( border ) < size;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV3(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilGridPos[], 
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_LOG_ERROR( logger,  "stencilGEMV on grid " << gridSizes[0] << " x " << gridSizes[1] << " x " << gridSizes[2] )

    const IndexType i_stripes[] = { 0, lb[0], gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, lb[1], gridSizes[1] - ub[1], gridSizes[1] };
    const IndexType k_stripes[] = { 0, lb[2], gridSizes[2] - ub[2], gridSizes[2] };

    /** The three dimensional grid is partitioned in each dimension according to the boundaries */

    for ( IndexType ib = 0; ib < 3; ib++ )
    {
        IndexType i0 = i_stripes[ib];
        IndexType i1 = i_stripes[ib + 1];

        for ( IndexType jb = 0; jb < 3; jb++ )
        {
            IndexType j0 = j_stripes[jb];
            IndexType j1 = j_stripes[jb + 1];

            for ( IndexType kb = 0; kb < 3; kb++ )
            {
                IndexType k0 = k_stripes[kb]; 
                IndexType k1 = k_stripes[kb + 1];

                SCAI_LOG_ERROR( logger, "loop over " << i0 << " - " << i1 << ", " << j0 << " - " << j1 << ", " << k0 << " - " << k1 )

                if ( ib == 1 && jb == 1 && kb == 1 )
                {
                    // this is the inner part, all stencil points are valid, no inner checks needed

                    #pragma omp parallel for

                    for ( IndexType i = i0; i < i1; i++ )
                    for ( IndexType j = j0; j < j1; j++ )
                    for ( IndexType k = k0; k < k1; k++ )
                    {
                        IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

                        for ( IndexType p = 0; p < nPoints; ++p )
                        {   
                            result[ gridPos ] += alpha * stencilVal[p] * x[ gridPos + stencilLinPos[p] ];
                        }
                    }
                }
                else
                {
                    // at least one stencil point neighbor might be invalid, check all


                    for ( IndexType i = i0; i < i1; i++ )
                    for ( IndexType j = j0; j < j1; j++ )
                    for ( IndexType k = k0; k < k1; k++ )
                    {
                        IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

                        for ( IndexType p = 0; p < nPoints; ++p )
                        {
                            if ( !isInner( i, stencilGridPos[ 3 * p ], gridSizes[0] ) )
                            {   
                                continue;
                            }
                            if ( !isInner( j, stencilGridPos[ 3 * p + 1 ], gridSizes[1] ) )
                            {   
                                continue;
                            }
                            if ( !isInner( k, stencilGridPos[ 3 * p + 2 ], gridSizes[2] ) )
                            {   
                                continue;
                            }
                            
                            result[ gridPos ] += alpha * stencilVal[p] * x[ gridPos + stencilLinPos[p] ];
                        }
                    }
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV( 
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims, 
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilGridPos[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    if ( nDims == 3 )
    {
        stencilGEMV3( result, alpha, x, gridSizes, lb, ub, gridDistances,
                      nPoints, stencilGridPos, stencilVal, stencilLinPos );
    }
    else
    {
        COMMON_THROWEXCEPTION( "stencilGEMV for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_DEBUG( logger,
                    "register StencilKernel OpenMP-routines for Host at kernel registry [" << flag 
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<StencilKernelTrait::stencilGEMV<ValueType> >( stencilGEMV, ctx, flag );

    // typename StencilKernelTrait::stencilGEMV<ValueType>::FuncType f = &stencilGEMV<ValueType>;
}

/* --------------------------------------------------------------------------- */

OpenMPStencilKernel::OpenMPStencilKernel()
{
    SCAI_LOG_INFO( logger, "register StencilKernel OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

OpenMPStencilKernel::~OpenMPStencilKernel()
{
    SCAI_LOG_INFO( logger, "unregister StencilKernel OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

OpenMPStencilKernel OpenMPStencilKernel::guard;

/* --------------------------------------------------------------------------- */

}

}
