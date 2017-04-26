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
#include <scai/tracing.hpp>

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

void OpenMPStencilKernel::stencilSizes1(
    IndexType sizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        IndexType gridPos = i * gridDistances[0];

        IndexType cnt = 0;

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, stencilNodes[ p ], gridSizes[0] ) )
            {   
                continue;
            }

            cnt++;

            sizes[ gridPos ] = cnt;
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilSizes2(
    IndexType sizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        for ( IndexType j = 0; j < gridSizes[1]; ++j )
        {
            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

            IndexType cnt = 0;

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if ( !isInner( i, stencilNodes[ 2 * p ], gridSizes[0] ) )
                {   
                    continue;
                }
                if ( !isInner( j, stencilNodes[ 2 * p + 1 ], gridSizes[1] ) )
                {   
                    continue;
                }
                cnt++;
            }

            sizes[ gridPos ] = cnt;
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilSizes3(
    IndexType sizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        for ( IndexType j = 0; j < gridSizes[1]; ++j )
        {
            for ( IndexType k = 0; k < gridSizes[2]; ++k )
            {
                IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

                IndexType cnt = 0;

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if ( !isInner( i, stencilNodes[ 3 * p ], gridSizes[0] ) )
                    {   
                        continue;
                    }
                    if ( !isInner( j, stencilNodes[ 3 * p + 1 ], gridSizes[1] ) )
                    {   
                        continue;
                    }
                    if ( !isInner( k, stencilNodes[ 3 * p + 2 ], gridSizes[2] ) )
                    {   
                        continue;
                    }
                    cnt++;
                }

                sizes[ gridPos ] = cnt;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilSizes(
    IndexType sizes[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    SCAI_REGION( "OpenMP.StencilSizes" )

    switch ( nDims ) 
    {
        case 1 : stencilSizes1( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        case 2 : stencilSizes2( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        case 3 : stencilSizes3( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilSizes for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencil2CSR1(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        IndexType gridPos = i * gridDistances[0];

        IndexType ia = csrIA[ gridPos ];

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, stencilNodes[ p ], gridSizes[0] ) )
            {   
                continue;
            }

            csrJA[ ia ] = gridPos + stencilLinPos[p];
            csrValues[ ia ] = stencilVal[p];

            ia++;
        }

        SCAI_ASSERT_EQ_ERROR( ia, csrIA[gridPos + 1], "serious mismatch, invalid csrIA array" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencil2CSR2(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        for ( IndexType j = 0; j < gridSizes[1]; ++j )
        {
            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

            IndexType ia = csrIA[ gridPos ];

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if ( !isInner( i, stencilNodes[ 2 * p ], gridSizes[0] ) )
                {   
                    continue;
                }
                if ( !isInner( j, stencilNodes[ 2 * p + 1 ], gridSizes[1] ) )
                {   
                    continue;
                }

                csrJA[ ia ] = gridPos + stencilLinPos[p];
                csrValues[ ia ] = stencilVal[p];

                ia++;
            }

            SCAI_ASSERT_EQ_ERROR( ia, csrIA[gridPos + 1], "serious mismatch, invalid csrIA array" )
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencil2CSR3(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        for ( IndexType j = 0; j < gridSizes[1]; ++j )
        {
            for ( IndexType k = 0; k < gridSizes[2]; ++k )
            {
                IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

                IndexType ia = csrIA[ gridPos ];

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if ( !isInner( i, stencilNodes[ 3 * p ], gridSizes[0] ) )
                    {   
                        continue;
                    }
                    if ( !isInner( j, stencilNodes[ 3 * p + 1 ], gridSizes[1] ) )
                    {   
                        continue;
                    }
                    if ( !isInner( k, stencilNodes[ 3 * p + 2 ], gridSizes[2] ) )
                    {   
                        continue;
                    }

                    csrJA[ ia ] = gridPos + stencilLinPos[p];
                    csrValues[ ia ] = stencilVal[p];

                    ia++;
                }

                SCAI_ASSERT_EQ_ERROR( ia, csrIA[gridPos + 1], "serious mismatch, invalid csrIA array" )
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencil2CSR(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_REGION( "OpenMP.Stencil2CSR" )

    switch ( nDims ) 
    {
        case 1 : stencil2CSR1( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 2 : stencil2CSR2( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 3 : stencil2CSR3( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        default: COMMON_THROWEXCEPTION( "stencil2CSR for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV1(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_LOG_INFO( logger,  "stencilGEMV on grid " << gridSizes[0] )

    const IndexType i_stripes[] = { 0, lb[0], gridSizes[0] - ub[0], gridSizes[0] };

    /** The two-dimensional grid is partitioned in each dimension according to the boundaries */

    for ( IndexType ib = 0; ib < 3; ib++ )
    {
        IndexType i0 = i_stripes[ib];
        IndexType i1 = i_stripes[ib + 1];

        SCAI_LOG_DEBUG( logger, "loop over " << i0 << " - " << i1 )

        if ( ib == 1 )
        {
            // this is the inner part, all stencil points are valid, no inner checks needed

            #pragma omp parallel for

            for ( IndexType i = i0; i < i1; i++ )
            {
                IndexType gridPos = i * gridDistances[0];

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
            {
                IndexType gridPos = i * gridDistances[0];

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if ( !isInner( i, stencilNodes[ p ], gridSizes[0] ) )
                    {   
                        continue;
                    }
                    
                    result[ gridPos ] += alpha * stencilVal[p] * x[ gridPos + stencilLinPos[p] ];
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV2(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_LOG_INFO( logger,  "stencilGEMV on grid " << gridSizes[0] << " x " << gridSizes[1] )

    const IndexType i_stripes[] = { 0, lb[0], gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, lb[1], gridSizes[1] - ub[1], gridSizes[1] };

    /** The two-dimensional grid is partitioned in each dimension according to the boundaries */

    for ( IndexType ib = 0; ib < 3; ib++ )
    {
        IndexType i0 = i_stripes[ib];
        IndexType i1 = i_stripes[ib + 1];

        for ( IndexType jb = 0; jb < 3; jb++ )
        {
            IndexType j0 = j_stripes[jb];
            IndexType j1 = j_stripes[jb + 1];

            SCAI_LOG_DEBUG( logger, "loop over " << i0 << " - " << i1 << ", " << j0 << " - " << j1 )

            if ( ib == 1 && jb == 1 )
            {
                // this is the inner part, all stencil points are valid, no inner checks needed

                #pragma omp parallel for

                for ( IndexType i = i0; i < i1; i++ )
                for ( IndexType j = j0; j < j1; j++ )
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

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
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( !isInner( i, stencilNodes[ 2 * p ], gridSizes[0] ) )
                        {   
                            continue;
                        }
                        if ( !isInner( j, stencilNodes[ 2 * p + 1 ], gridSizes[1] ) )
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
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_LOG_INFO( logger,  "stencilGEMV on grid " << gridSizes[0] << " x " << gridSizes[1] << " x " << gridSizes[2] )

    const IndexType nStripes = 3;

    const IndexType i_stripes[] = { 0, lb[0], gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, lb[1], gridSizes[1] - ub[1], gridSizes[1] };
    const IndexType k_stripes[] = { 0, lb[2], gridSizes[2] - ub[2], gridSizes[2] };

    /** The three dimensional grid is partitioned in each dimension according to the boundaries */

    for ( IndexType ib = 0; ib < nStripes; ib++ )
    {
        IndexType i0 = i_stripes[ib];
        IndexType i1 = i_stripes[ib + 1];

        for ( IndexType jb = 0; jb < nStripes; jb++ )
        {
            IndexType j0 = j_stripes[jb];
            IndexType j1 = j_stripes[jb + 1];

            for ( IndexType kb = 0; kb < nStripes; kb++ )
            {
                IndexType k0 = k_stripes[kb]; 
                IndexType k1 = k_stripes[kb + 1];

                SCAI_LOG_DEBUG( logger, "loop over " << i0 << " - " << i1 << ", " << j0 << " - " << j1 << ", " << k0 << " - " << k1 )

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
                            if ( !isInner( i, stencilNodes[ 3 * p ], gridSizes[0] ) )
                            {   
                                continue;
                            }
                            if ( !isInner( j, stencilNodes[ 3 * p + 1 ], gridSizes[1] ) )
                            {   
                                continue;
                            }
                            if ( !isInner( k, stencilNodes[ 3 * p + 2 ], gridSizes[2] ) )
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
void OpenMPStencilKernel::stencilGEMV4(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_LOG_INFO( logger,  "stencilGEMV on grid " << gridSizes[0] << " x " << gridSizes[1] 
                                           << " x " << gridSizes[2] << " x " << gridSizes[3] )

    const IndexType nStripes = 4;

    const IndexType i_stripes[] = { 0, lb[0], gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, lb[1], gridSizes[1] - ub[1], gridSizes[1] };
    const IndexType k_stripes[] = { 0, lb[2], gridSizes[2] - ub[2], gridSizes[2] };
    const IndexType m_stripes[] = { 0, lb[3], gridSizes[3] - ub[3], gridSizes[3] };

    /** The four-dimensional grid is partitioned in each dimension according to the boundaries */

    for ( IndexType ib = 0; ib < nStripes; ib++ )
    {
        IndexType i0 = i_stripes[ib];
        IndexType i1 = i_stripes[ib + 1];

        for ( IndexType jb = 0; jb < nStripes; jb++ )
        {
            IndexType j0 = j_stripes[jb];
            IndexType j1 = j_stripes[jb + 1];

            for ( IndexType kb = 0; kb < nStripes; kb++ )
            {
                IndexType k0 = k_stripes[kb]; 
                IndexType k1 = k_stripes[kb + 1];

                for ( IndexType mb = 0; mb < nStripes; mb++ )
                {
                    IndexType m0 = m_stripes[mb]; 
                    IndexType m1 = m_stripes[mb + 1];

                    SCAI_LOG_DEBUG( logger, "loop over " << i0 << " - " << i1 << ", " << j0 << " - " << j1 
                                                 << ", " << k0 << " - " << k1 << ", " << m0 << " - " << m1 )

                    if ( ib == 1 && jb == 1 && kb == 1 && mb == 1 )
                    {
                        // this is the inner part, all stencil points are valid, no inner checks needed

                        #pragma omp parallel for

                        for ( IndexType i = i0; i < i1; i++ )
                        for ( IndexType j = j0; j < j1; j++ )
                        for ( IndexType k = k0; k < k1; k++ )
                        for ( IndexType m = m0; m < m1; m++ )
                        {
                            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

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
                        for ( IndexType m = m0; m < m1; m++ )
                        {
                            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];
    
                            for ( IndexType p = 0; p < nPoints; ++p )
                            {
                                if ( !isInner( i, stencilNodes[ 4 * p ], gridSizes[0] ) )
                                {   
                                    continue;
                                }
                                if ( !isInner( j, stencilNodes[ 4 * p + 1 ], gridSizes[1] ) )
                                {   
                                    continue;
                                }
                                if ( !isInner( k, stencilNodes[ 4 * p + 2 ], gridSizes[2] ) )
                                {   
                                    continue;
                                }
                                if ( !isInner( m, stencilNodes[ 4 * p + 3 ], gridSizes[3] ) )
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
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV" )

    switch ( nDims ) 
    {
        case 1 : stencilGEMV1( result, alpha, x, gridSizes, lb, ub, gridDistances,
                               nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 2 : stencilGEMV2( result, alpha, x, gridSizes, lb, ub, gridDistances,
                               nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 3 : stencilGEMV3( result, alpha, x, gridSizes, lb, ub, gridDistances,
                               nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 4 : stencilGEMV4( result, alpha, x, gridSizes, lb, ub, gridDistances,
                               nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMV for nDims = " << nDims << " not supported yet" )
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

    KernelRegistry::set<StencilKernelTrait::stencil2CSR<ValueType> >( stencil2CSR, ctx, flag );
    KernelRegistry::set<StencilKernelTrait::stencilGEMV<ValueType> >( stencilGEMV, ctx, flag );
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
