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

#include <scai/common/Grid.hpp>
#include <scai/sparsekernel/openmp/OpenMPStencilKernel.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>

namespace scai
{
namespace sparsekernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPStencilKernel::logger, "OpenMP.StencilKernel" )

/* --------------------------------------------------------------------------- */

/** This routine checks whether a point pos is an inner point 
 *
 *  @param[in] pos is the position, in range 0 .. size-1
 *  @param[in] offset specifies the relative offset to pos, might be positive or negative
 *  @param[in] size is the size of the range, only needed for offset > 0
 */
static inline bool isInner( const IndexType pos, const int offset, const IndexType size )
{
    if ( offset == 0 )
    {
        return true;
    }
    if ( offset < 0 )
    {
        return pos >= static_cast<IndexType>( -offset );
    }

    return pos + static_cast<IndexType>( offset ) < size;
}

/** Inline predicate to check if a stencil point is still in the grid */

static inline bool isInner2( const IndexType pos0, const IndexType pos1, 
                             const int offset[2], const IndexType size[2] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    return true;
}

/** Inline predicate to check if a stencil point is still in the grid */

static inline bool isInner3( const IndexType pos0, const IndexType pos1, const IndexType pos2, 
                             const int offset[3], const IndexType size[3] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    if ( !isInner( pos2, offset[2], size[2] ) ) return false;
    return true;
}

/** Inline predicate to check if a 4-dimensional stencil point is still in the grid */

static inline bool isInner4( const IndexType pos0, const IndexType pos1, 
                             const IndexType pos2, const IndexType pos3, 
                             const int offset[4], const IndexType size[4] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    if ( !isInner( pos2, offset[2], size[2] ) ) return false;
    if ( !isInner( pos3, offset[3], size[3] ) ) return false;
    return true;
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
                if ( isInner2( i, j, &stencilNodes[ 2 * p ], gridSizes ) )
                {   
                    cnt++;
                }
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
                    if ( isInner3( i, j, k, &stencilNodes[ 3 * p ], gridSizes ) )
                    {   
                        cnt++;
                    }
                }

                sizes[ gridPos ] = cnt;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilLocalSizes4(
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
                for ( IndexType m = 0; m < gridSizes[3]; ++m )
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

                    IndexType cnt = 0;

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( isInner4( i, j, k, m, &stencilNodes[ 4 * p ], gridSizes ) )
                        { 
                            cnt++;
                        }
                    }
    
                    sizes[ gridPos ] = cnt;
                }
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

        case 4 : stencilLocalSizes4( sizes, gridSizes, gridDistances, nPoints, stencilNodes );
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
    const int stencilOffset[] )
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

            csrJA[ ia ] = gridPos + stencilOffset[p];
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
    const int stencilOffset[] )
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
                if ( isInner2( i, j, &stencilNodes[2 * p], gridSizes ) )
                {   
                    csrJA[ ia ] = gridPos + stencilOffset[p];
                    csrValues[ ia ] = stencilVal[p];
                    ia++;
                }
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
    const int stencilOffset[] )
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
                    if ( isInner3( i, j, k, &stencilNodes[3 * p], gridSizes ) )
                    {   
                        csrJA[ ia ] = gridPos + stencilOffset[p];
                        csrValues[ ia ] = stencilVal[p];
                        ia++;
                    }
                }

                SCAI_ASSERT_EQ_ERROR( ia, csrIA[gridPos + 1], "serious mismatch, invalid csrIA array" )
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilLocalCSR4(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        for ( IndexType j = 0; j < gridSizes[1]; ++j )
        {
            for ( IndexType k = 0; k < gridSizes[2]; ++k )
            {
                for ( IndexType m = 0; m < gridSizes[3]; ++m )
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

                    IndexType ia = csrIA[ gridPos ];

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( isInner4( i, j, k, m, &stencilNodes[ 4 * p ], gridSizes ) )
                        {   
                            csrJA[ ia ] = gridPos + stencilOffset[p];
                            csrValues[ ia ] = stencilVal[p];
    
                            ia++;
                        }
                    }

                    SCAI_ASSERT_EQ_ERROR( ia, csrIA[gridPos + 1], "serious mismatch, invalid csrIA array" )
                }
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
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.LocalCSR" )

    switch ( nDims ) 
    {
        case 1 : stencilLocalCSR1( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilLocalCSR2( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilLocalCSR3( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilLocalCSR4( csrJA, csrValues, csrIA, gridSizes, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilLocalCSR for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes1(
    IndexType sizes[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    // traverse over all points of the local grid

    #pragma omp parallel for
    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType iGlobal = localLB[0] + iLocal;

        IndexType localIndex = iLocal * localGridDistances[0];

        IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

        // check for each stencil point if it is not local

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( isInner( iLocal, stencilNodes[p], localGridSizes[0] ) )
            {
                // stencil point is local, no halo entry
                continue;
            }
            if ( isInner( iGlobal, stencilNodes[p], globalGridSizes[0] ) )
            {
                // stencil point is not local, but globally available, so count it
                cnt ++;
            }
        }

        sizes[ localIndex ] = cnt;
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

    #pragma omp parallel for
    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType iGlobal = localLB[0] + iLocal;

        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            IndexType localIndex = iLocal * localGridDistances[0] + jLocal * localGridDistances[1];

            IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

            // check for each stencil point if it is not local

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if ( isInner2( iLocal, jLocal, &stencilNodes[2 * p], localGridSizes ) )
                {
                    // stencil point is local, no halo entry
                    continue;
                }
                if ( isInner2( iGlobal, jGlobal, &stencilNodes[2 * p], globalGridSizes ) )
                {
                    // stencil point is not local, but globally available, so count it
                    cnt ++;
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

    #pragma omp parallel for

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType iGlobal = localLB[0] + iLocal;

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
                    if ( isInner3( iLocal, jLocal, kLocal, &stencilNodes[ 3 * p ], localGridSizes ) )
                    {
                        // stencil point is local, no halo entry
                        continue;
                    }
                    if ( isInner3( iGlobal, jGlobal, kGlobal, &stencilNodes[3 * p], globalGridSizes ) )
                    {
                        // stencil point is not local, but globally available, so count it
                        cnt ++;
                    }
                }

                sizes[ localIndex ] = cnt;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes4(
    IndexType sizes[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    // traverse over all points of the local 4-dimensional grid

    #pragma omp parallel for

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType iGlobal = localLB[0] + iLocal;

        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            for ( IndexType kLocal = 0, kGlobal = localLB[2]; kLocal < localGridSizes[2]; ++kLocal, ++kGlobal )
            {
                for ( IndexType mLocal = 0, mGlobal = localLB[3]; mLocal < localGridSizes[3]; ++mLocal, ++mGlobal )
                {
                    IndexType localIndex =   iLocal * localGridDistances[0] + jLocal * localGridDistances[1] 
                                           + kLocal * localGridDistances[2] + mLocal * localGridDistances[3];
    
                    IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

                    // check for each stencil point if it is not local

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( isInner4( iLocal, jLocal, kLocal, mLocal, &stencilNodes[4 * p], localGridSizes ) )
                        {
                            continue;   // local stencil neighbor, does not count for halo
                        }
                        if ( isInner4( iGlobal, jGlobal, kGlobal, mGlobal, &stencilNodes[4 * p], globalGridSizes ) )
                        {
                            cnt ++;     // stencil point is not local, but globally available, so count it
                        }
                    }

                    sizes[ localIndex ] = cnt;
                }
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
        case 1 : stencilHaloSizes1( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        case 2 : stencilHaloSizes2( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        case 3 : stencilHaloSizes3( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        case 4 : stencilHaloSizes4( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, nPoints, stencilNodes );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilHaloSizes for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilHaloCSR1(
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
    const int stencilOffset[] )
{
    // traverse over all points of the local grid

    #pragma omp parallel for

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType localIndex  = iLocal * localGridDistances[0];

        IndexType offset = csrIA[ localIndex ];

        if ( csrIA[ localIndex + 1 ] == offset )
        {
            continue;   // this grid point has all its relevant neighbors locally, no halo
        }

        IndexType iGlobal = localLB[0] + iLocal;
        IndexType globalIndex = iGlobal * globalGridDistances[0];

        // check for each stencil point if it is not local

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( isInner( iLocal, stencilNodes[p], localGridSizes[0] ) )
            {
                continue;      // stencil point is local, no halo entry
            }

            // stencil point must still be globally available to be counted

            if ( isInner( iGlobal, stencilNodes[p], globalGridSizes[0] ) )
            {
                csrJA[offset] = globalIndex + stencilOffset[p];
                csrValues[offset] = stencilVal[p];
                ++offset;
            }
        }

        SCAI_ASSERT_EQ_ERROR( offset, csrIA[ localIndex + 1 ], "serious mismatch" );
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
    const int stencilOffset[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        const IndexType iGlobal = localLB[0] + iLocal;

        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            IndexType localIndex  = iLocal * localGridDistances[0] + jLocal * localGridDistances[1];

            IndexType offset = csrIA[ localIndex ];

            if ( csrIA[ localIndex + 1 ] == offset )
            {
                continue;   // this grid point has all its relevant neighbors locally, no halo
            }

            IndexType globalIndex = iGlobal * globalGridDistances[0] + jGlobal * globalGridDistances[1];

            // check for each stencil point if it is not local

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if ( isInner2( iLocal, jLocal, &stencilNodes[2 * p], localGridSizes ) )
                {
                    continue;      // stencil point is local, no halo entry
                }

                // stencil point must still be globally available to be counted

                if ( isInner2( iGlobal, jGlobal, &stencilNodes[2 * p], globalGridSizes ) )
                {
                    csrJA[offset] = globalIndex + stencilOffset[p];
                    csrValues[offset] = stencilVal[p];
                    ++offset;
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
    const int stencilOffset[] )
{
    // traverse over all points of the local grid

    #pragma omp parallel for

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        const IndexType iGlobal = localLB[0] + iLocal;

        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            for ( IndexType kLocal = 0, kGlobal = localLB[2]; kLocal < localGridSizes[2]; ++kLocal, ++kGlobal )
            {
                IndexType localIndex =   iLocal * localGridDistances[0] 
                                       + jLocal * localGridDistances[1] 
                                       + kLocal * localGridDistances[2];

                IndexType offset = csrIA[ localIndex ];

                if ( csrIA[ localIndex + 1 ] == offset )
                {
                    continue;   // this grid point has all its relevant neighbors locally, no halo
                }

                IndexType globalIndex =    iGlobal * globalGridDistances[0] 
                                         + jGlobal * globalGridDistances[1]
                                         + kGlobal * globalGridDistances[2];

                // check for each stencil point if it is not local

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if ( isInner3( iLocal, jLocal, kLocal, &stencilNodes[ 3 * p ], localGridSizes ) )
                    {
                        // stencil point is local, no halo entry
                        continue;
                    }

                    // stencil point must still be globally available to be counted

                    if ( isInner3( iGlobal, jGlobal, kGlobal, &stencilNodes[3 * p], globalGridSizes ) )
                    {
                        csrJA[offset] = globalIndex + stencilOffset[p];
                        csrValues[offset] = stencilVal[p];
                        ++offset;
                    }
                }
     
                SCAI_ASSERT_EQ_ERROR( offset, csrIA[ localIndex + 1 ], "serious mismatch" );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilHaloCSR4(
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
    const int stencilOffset[] )
{
    // traverse over all points of the local grid

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        const IndexType iGlobal = localLB[0] + iLocal;

        for ( IndexType jLocal = 0, jGlobal = localLB[1]; jLocal < localGridSizes[1]; ++jLocal, ++jGlobal )
        {
            for ( IndexType kLocal = 0, kGlobal = localLB[2]; kLocal < localGridSizes[2]; ++kLocal, ++kGlobal )
            {
                for ( IndexType mLocal = 0, mGlobal = localLB[3]; mLocal < localGridSizes[3]; ++mLocal, ++mGlobal )
                {
                    IndexType localIndex =   iLocal * localGridDistances[0] 
                                           + jLocal * localGridDistances[1] 
                                           + kLocal * localGridDistances[2]
                                           + mLocal * localGridDistances[3];

                    IndexType offset = csrIA[ localIndex ];

                    if ( csrIA[ localIndex + 1 ] == offset )
                    {
                        continue;   // this grid point has all its relevant neighbors locally, no halo
                    }

                    IndexType globalIndex =    iGlobal * globalGridDistances[0] 
                                             + jGlobal * globalGridDistances[1]
                                             + kGlobal * globalGridDistances[2]
                                             + mGlobal * globalGridDistances[3];

                    // check for each stencil point if it is not local

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( isInner4( iGlobal, jGlobal, kGlobal, mGlobal, &stencilNodes[4 * p], globalGridSizes ) )
                        {
                            // stencil point is in global grid
    
                            if ( !isInner4( iLocal, jLocal, kLocal, mLocal, &stencilNodes[4 * p], localGridSizes ) )
                            {
                                // but not locally, so it belongs to halo, store global index

                                csrJA[offset] = globalIndex + stencilOffset[p];
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
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.HaloCSR" )

    switch ( nDims ) 
    {
        case 1 : stencilHaloCSR1( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilHaloCSR2( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilHaloCSR3( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilHaloCSR4( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilHaloCSR for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV1Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV1Inner" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];

    SCAI_LOG_INFO( logger,  "stencilGEMV1Inner on " << i0 << " - " << i1 )

    #pragma omp parallel for
#if defined( __INTEL_COMPILER )
    #pragma unroll_and_jam( 4 )
#endif
    for ( IndexType i = i0; i < i1; i++ )
    {
        IndexType gridPos = i * gridDistances[0];

        ValueType v = 0;

        for ( IndexType p = 0; p < nPoints; ++p )
        {   
            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
        }

        result[ gridPos] += alpha * v;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV2Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV2Inner" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];

    SCAI_LOG_INFO( logger,  "stencilGEMV2Inner on " << i0 << " - " << i1 << " x " << j0 << " - " << j1 )

    #pragma omp parallel for
    for ( IndexType i = i0; i < i1; i++ )
    {
#if defined( __INTEL_COMPILER )
        #pragma unroll_and_jam( 4 )
#endif
        for ( IndexType j = j0; j < j1; j++ )
        {
            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

            ValueType v = 0;

            for ( IndexType p = 0; p < nPoints; ++p )
            {   
                v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
            }

            result[ gridPos] += alpha * v;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV3Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV3Inner" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];

    SCAI_LOG_INFO( logger,  "stencilGEMV3Inner on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 
                             << " x " << k0 << " - " << k1  )

    #pragma omp parallel for
    for ( IndexType i = i0; i < i1; i++ )
    {
        for ( IndexType j = j0; j < j1; j++ )
        {
#if defined( __INTEL_COMPILER )
            #pragma unroll_and_jam( 4 )
#endif
            for ( IndexType k = k0; k < k1; k++ )
            {
                IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];
    
                ValueType v = 0;
    
                for ( IndexType p = 0; p < nPoints; ++p )
                {   
                    v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
                }
    
                result[ gridPos] += alpha * v;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV4Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV4Inner" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];
    const IndexType m0 = gridBounds[6];
    const IndexType m1 = gridBounds[7];

    SCAI_LOG_INFO( logger,  "stencilGEMV4Inner ( " << nPoints << " points ) on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 
                             << " x " << k0 << " - " << k1  
                             << " x " << m0 << " - " << m1  )

    #pragma omp parallel for
    for ( IndexType i = i0; i < i1; i++ )
    {
        for ( IndexType j = j0; j < j1; j++ )
        {
            for ( IndexType k = k0; k < k1; k++ )
            {
#if defined( __INTEL_COMPILER )
                #pragma unroll_and_jam( 4 )
#endif
                for ( IndexType m = m0; m < m1; m++ )
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];
    
                    ValueType v = 0;
        
                    for ( IndexType p = 0; p < nPoints; ++p )
                    {   
                        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
                    }
    
                    result[ gridPos] += alpha * v;
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMVInner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    switch ( nDims ) 
    {
        case 1 : stencilGEMV1Inner( result, alpha, x, gridBounds, gridDistances,
                                    nPoints, stencilVal, stencilOffset );
                 break;

        case 2 : stencilGEMV2Inner( result, alpha, x, gridBounds, gridDistances,
                                    nPoints, stencilVal, stencilOffset );
                 break;

        case 3 : stencilGEMV3Inner( result, alpha, x, gridBounds, gridDistances,
                                    nPoints, stencilVal, stencilOffset );
                 break;

        case 4 : stencilGEMV4Inner( result, alpha, x, gridBounds, gridDistances,
                                    nPoints, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMVInner for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV1Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV1Border" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];

    SCAI_LOG_INFO( logger,  "stencilGEMV1Border on " << i0 << " - " << i1 )

    #pragma omp parallel for if ( i1 - i0 > 7 )
    for ( IndexType i = i0; i < i1; i++ )
    {
        IndexType gridPos = i * gridDistances[0];

        ValueType v = 0;

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( isInner( i, stencilNodes[ p ], gridSizes [ 0 ] ) )
            {
                v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
            }
        }

        result[ gridPos] += alpha * v;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV2Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV2Border" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];

    SCAI_LOG_INFO( logger,  "stencilGEMV2Border on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 )

    bool outer_parallel = ( i1 - i0 ) > 7;

    #pragma omp parallel for if ( outer_parallel )
    for ( IndexType i = i0; i < i1; i++ )
    {
        #pragma omp parallel for if ( !outer_parallel )
        for ( IndexType j = j0; j < j1; j++ )
        {
            IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

            ValueType v = 0;

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                if ( isInner2( i, j, &stencilNodes[ 2 * p ], gridSizes ) )
                {
                    v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
                }
            }

            result[ gridPos] += alpha * v;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV3Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV3Border" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];

    SCAI_LOG_INFO( logger,  "stencilGEMV3Border on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 
                             << " x " << k0 << " - " << k1  )

    bool outer_parallel = ( i1 - i0 ) > 7 ;

    #pragma omp parallel for if ( outer_parallel )
    for ( IndexType i = i0; i < i1; i++ )
    {
        #pragma omp parallel for if ( !outer_parallel )
        for ( IndexType j = j0; j < j1; j++ )
        {
            for ( IndexType k = k0; k < k1; k++ )
            {
                IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

                ValueType v = 0;

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    if ( isInner3( i, j, k, &stencilNodes[ 3 * p ], gridSizes ) )
                    {
                        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
                    }
                }

                result[ gridPos] += alpha * v;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV4Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV4Border" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];
    const IndexType m0 = gridBounds[6];
    const IndexType m1 = gridBounds[7];

    SCAI_LOG_INFO( logger,  "stencilGEMV4Border on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 
                             << " x " << k0 << " - " << k1  
                             << " x " << m0 << " - " << m1  )

    // take decision which of the two nested loops to take for OpenMP parallelization

    bool outer_parallel = ( i1 - i0 ) > 7 ;

    #pragma omp parallel for if ( outer_parallel )
    for ( IndexType i = i0; i < i1; i++ )
    {
        #pragma omp parallel for if ( !outer_parallel )
        for ( IndexType j = j0; j < j1; j++ )
        {
            for ( IndexType k = k0; k < k1; k++ )
            {
                for ( IndexType m = m0; m < m1; m++ )
                {
                    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

                    ValueType v = 0;

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        if ( isInner4( i, j, k, m, &stencilNodes[ 4 * p ], gridSizes ) )
                        {
                            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
                        }
                    }

                    result[ gridPos] += alpha * v;
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMVBorder(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    switch ( nDims ) 
    {
        case 1 : stencilGEMV1Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilGEMV2Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilGEMV3Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilGEMV4Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMVBorder for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMVCaller(
    IndexType gridBounds[],
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType width[],
    const IndexType gridDistances[],
    const IndexType currentDim,
    const IndexType nPoints,
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( currentDim, nDims, "illegal current dimension" )

    if ( gridSizes[currentDim] <= width[2 * currentDim] + width[2 * currentDim + 1] )
    {
        // no inner part in current dimension, call boundary routine, afterwards all is done

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
 
        return;
    }

    // if left boundary, handle it

    if ( width[2 * currentDim] > 0 )
    {
        gridBounds[2 * currentDim] = 0;
        gridBounds[2 * currentDim + 1] = width[2 * currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // set inner boundaries,

    gridBounds[2 * currentDim] = width[2 * currentDim];
    gridBounds[2 * currentDim + 1] = gridSizes[currentDim] - width[2 * currentDim + 1];

    if ( currentDim + 1 == nDims )
    {
        // all boundaries are now inner ones, we can call the routine for the inner points

        stencilGEMVInner( result, alpha, x, nDims, gridBounds, gridDistances, nPoints, stencilVal, stencilOffset );
    }
    else
    { 
        // recursive call to set grid bounds for the next dimension 

        stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, width, gridDistances, currentDim + 1, 
                           nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // if right boundary, travere it

    if ( width[2 * currentDim + 1] > 0 )
    {
        gridBounds[2 * currentDim] = gridSizes[currentDim] - width[2 * currentDim + 1];
        gridBounds[2 * currentDim + 1] = gridSizes[currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // reset the boundaries

    gridBounds[2 * currentDim] = 0;
    gridBounds[2 * currentDim + 1] = gridSizes[currentDim];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::stencilGEMV( 
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims, 
    const IndexType gridSizes[],
    const IndexType width[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    // prepare array with gridBounds for the recursive traverser routine

    IndexType gridBounds[ 2 * SCAI_GRID_MAX_DIMENSION ];

    for ( IndexType i = 0; i < nDims; ++i )
    {
        gridBounds[2 * i] = 0;
        gridBounds[2 * i + 1] = gridSizes[i];
    }

    IndexType currentDim = 0;

    stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, width, gridDistances, currentDim,
                       nPoints, stencilNodes, stencilVal, stencilOffset );
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry; 

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_DEBUG( logger, "register StencilKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<StencilKernelTrait::stencilLocalSizes>( stencilLocalSizes, ctx, flag );
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

    Registrator::registerKernels( flag );

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

OpenMPStencilKernel::~OpenMPStencilKernel()
{
    SCAI_LOG_INFO( logger, "unregister StencilKernel OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::registerKernels( flag );

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

OpenMPStencilKernel OpenMPStencilKernel::guard;

/* --------------------------------------------------------------------------- */

}

}
