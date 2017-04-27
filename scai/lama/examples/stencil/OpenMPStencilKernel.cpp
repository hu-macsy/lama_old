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

void OpenMPStencilKernel::stencilLocalSizes1(
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

void OpenMPStencilKernel::stencilLocalSizes2(
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

void OpenMPStencilKernel::stencilLocalSizes3(
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

void OpenMPStencilKernel::stencilLocalSizes(
    IndexType sizes[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    SCAI_REGION( "OpenMP.Stencil.LocalSizes" )

    switch ( nDims ) 
    {
        case 1 : stencilLocalSizes1( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        case 2 : stencilLocalSizes2( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        case 3 : stencilLocalSizes3( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilLocalSizes for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilLocalCSR1(
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
void OpenMPStencilKernel::stencilLocalCSR2(
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
void OpenMPStencilKernel::stencilLocalCSR3(
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
void OpenMPStencilKernel::stencilLocalCSR(
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
    SCAI_REGION( "OpenMP.Stencil.LocalCSR" )

    switch ( nDims ) 
    {
        case 1 : stencilLocalCSR1( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 2 : stencilLocalCSR2( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 3 : stencilLocalCSR3( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilLocalCSR for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes2(
    IndexType sizes[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0, iGlobal = localLB[0]; iLocal < localGridSizes[0]; ++iLocal, ++iGlobal )
    {
        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            IndexType localIndex = iLocal * localGridDistances[0] + jLocal * localGridDistances[1];

            IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

            // check for each stencil point if it is not local

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if (    isInner( iGlobal, stencilNodes[2 * p], globalGridSizes[0] ) 
                     && isInner( jGlobal, stencilNodes[2 * p + 1], globalGridSizes[1] ) )
                {
                    // stencil point is in global grid

                    if (    isInner( iLocal, stencilNodes[ 2 * p ], localGridSizes[0] )
                         && isInner( jLocal, stencilNodes[ 2 * p + 1 ], localGridSizes[1] ) )
                    {
                        // stencil point is in local grid, do not count
                    }
                    else
                    {
                        cnt ++;
                    }
                }
            }

            sizes[ localIndex ] = cnt;
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes3(
    IndexType sizes[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0, iGlobal = localLB[0]; iLocal < localGridSizes[0]; ++iLocal, ++iGlobal )
    {
        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            for ( IndexType kLocal = 0, kGlobal = localLB[2]; kLocal < localGridSizes[2]; ++kLocal, ++kGlobal )
            {
                IndexType localIndex = iLocal * localGridDistances[0] + jLocal * localGridDistances[1] 
                                       + kLocal * localGridDistances[2];

                IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

                // check for each stencil point if it is not local

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if (    isInner( iGlobal, stencilNodes[3 * p], globalGridSizes[0] ) 
                         && isInner( jGlobal, stencilNodes[3 * p + 1], globalGridSizes[1] ) 
                         && isInner( kGlobal, stencilNodes[3 * p + 2], globalGridSizes[2] ) )
                    {
                        // stencil point is in global grid
    
                        if (    isInner( iLocal, stencilNodes[ 3 * p ], localGridSizes[0] )
                             && isInner( jLocal, stencilNodes[ 3 * p + 1 ], localGridSizes[1] ) 
                             && isInner( kLocal, stencilNodes[ 3 * p + 2 ], localGridSizes[2] ) )
                        {
                            // stencil point is in local grid, do not count
                        }
                        else
                        {
                            cnt ++;
                        }
                    }
                }

                sizes[ localIndex ] = cnt;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes(
    IndexType sizes[],
    const IndexType nDims,
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    SCAI_REGION( "OpenMP.Stencil.HaloSizes" )

    switch ( nDims ) 
    {
        case 2 : stencilHaloSizes2( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        case 3 : stencilHaloSizes3( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilLocalSizes for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilHaloCSR2(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType globalGridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0, iGlobal = localLB[0]; iLocal < localGridSizes[0]; ++iLocal, ++iGlobal )
    {
        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            IndexType localIndex  = iLocal * localGridDistances[0] + jLocal * localGridDistances[1];
            IndexType globalIndex = iGlobal * globalGridDistances[0] + jGlobal * globalGridDistances[1];

            IndexType offset = csrIA[ localIndex ];

            // check for each stencil point if it is not local

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if (    isInner( iGlobal, stencilNodes[2 * p], globalGridSizes[0] ) 
                     && isInner( jGlobal, stencilNodes[2 * p + 1], globalGridSizes[1] ) )
                {
                    // stencil point is in global grid

                    if (    isInner( iLocal, stencilNodes[ 2 * p ], localGridSizes[0] )
                         && isInner( jLocal, stencilNodes[ 2 * p + 1 ], localGridSizes[1] ) )
                    {
                        // stencil point is in local grid, do not count
                    }
                    else
                    {
                        csrJA[offset] = globalIndex + stencilLinPos[p];
                        csrValues[offset] = stencilVal[p];
                        ++offset;
                    }
                }
            }
 
            SCAI_ASSERT_EQ_ERROR( offset, csrIA[ localIndex + 1 ], "serious mismatch" );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilHaloCSR3(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType globalGridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0, iGlobal = localLB[0]; iLocal < localGridSizes[0]; ++iLocal, ++iGlobal )
    {
        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            for ( IndexType kLocal = 0, kGlobal = localLB[2]; kLocal < localGridSizes[2]; ++kLocal, ++kGlobal )
            {
                IndexType localIndex  = iLocal * localGridDistances[0] + jLocal * localGridDistances[1] 
                                         + kLocal * localGridDistances[2];
                IndexType globalIndex = iGlobal * globalGridDistances[0] + jGlobal * globalGridDistances[1]
                                         + kGlobal * globalGridDistances[2];

                IndexType offset = csrIA[ localIndex ];

                // check for each stencil point if it is not local

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if (    isInner( iGlobal, stencilNodes[3 * p], globalGridSizes[0] ) 
                         && isInner( jGlobal, stencilNodes[3 * p + 1], globalGridSizes[1] ) 
                         && isInner( kGlobal, stencilNodes[3 * p + 2], globalGridSizes[2] ) )
                    {
                        // stencil point is in global grid
    
                        if (    isInner( iLocal, stencilNodes[ 3 * p ], localGridSizes[0] )
                             && isInner( jLocal, stencilNodes[ 3 * p + 1 ], localGridSizes[1] ) 
                             && isInner( kLocal, stencilNodes[ 3 * p + 2 ], localGridSizes[2] ) )
                        {
                            // stencil point is in local grid, do not count
                        }
                        else
                        {
                            csrJA[offset] = globalIndex + stencilLinPos[p];
                            csrValues[offset] = stencilVal[p];
                            ++offset;
                        }
                    }
                }
     
                SCAI_ASSERT_EQ_ERROR( offset, csrIA[ localIndex + 1 ], "serious mismatch" );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilHaloCSR(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType nDims,
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType globalGridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilLinPos[] )
{
    SCAI_REGION( "OpenMP.Stencil.HaloCSR" )

    switch ( nDims ) 
    {
        case 2 : stencilHaloCSR2( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        case 3 : stencilHaloCSR3( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilLinPos );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilHaloCSR for nDims = " << nDims << " not supported yet" )
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

    bool small0 = gridSizes[0] <= lb[0] + ub[0];  // no inner part

    const IndexType i_stripes[] = { 0, small0 ? 0 : lb[0], small0 ? 0 : gridSizes[0] - ub[0], gridSizes[0] };

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

    bool small0 = gridSizes[0] <= lb[0] + ub[0];  // no inner part
    bool small1 = gridSizes[1] <= lb[1] + ub[1];  // no inner part

    const IndexType i_stripes[] = { 0, small0 ? 0 : lb[0], small0 ? 0 : gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, small1 ? 0 : lb[1], small1 ? 0 : gridSizes[1] - ub[1], gridSizes[1] };

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

    bool small0 = gridSizes[0] <= lb[0] + ub[0];  // no inner part
    bool small1 = gridSizes[1] <= lb[1] + ub[1];  // no inner part
    bool small2 = gridSizes[2] <= lb[2] + ub[2];  // no inner part

    const IndexType i_stripes[] = { 0, small0 ? 0 : lb[0], small0 ? 0 : gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, small1 ? 0 : lb[1], small1 ? 0 : gridSizes[1] - ub[1], gridSizes[1] };
    const IndexType k_stripes[] = { 0, small2 ? 0 : lb[2], small2 ? 0 : gridSizes[2] - ub[2], gridSizes[2] };

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

    bool small0 = gridSizes[0] <= lb[0] + ub[0];  // no inner part
    bool small1 = gridSizes[1] <= lb[1] + ub[1];  // no inner part
    bool small2 = gridSizes[2] <= lb[2] + ub[2];  // no inner part
    bool small3 = gridSizes[3] <= lb[3] + ub[3];  // no inner part

    const IndexType i_stripes[] = { 0, small0 ? 0 : lb[0], small0 ? 0 : gridSizes[0] - ub[0], gridSizes[0] };
    const IndexType j_stripes[] = { 0, small1 ? 0 : lb[1], small1 ? 0 : gridSizes[1] - ub[1], gridSizes[1] };
    const IndexType k_stripes[] = { 0, small2 ? 0 : lb[2], small2 ? 0 : gridSizes[2] - ub[2], gridSizes[2] };
    const IndexType m_stripes[] = { 0, small3 ? 0 : lb[3], small3 ? 0 : gridSizes[3] - ub[3], gridSizes[3] };

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

    KernelRegistry::set<StencilKernelTrait::stencilLocalCSR<ValueType> >( stencilLocalCSR, ctx, flag );
    KernelRegistry::set<StencilKernelTrait::stencilHaloCSR<ValueType> >( stencilHaloCSR, ctx, flag );
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
