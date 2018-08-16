/**
 * @file OpenMPStencilKernel.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief OpenMP Implementations on host for stencil kernels.
 * @author Thomas Brandes
 * @date 28.04.2017
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

void OpenMPStencilKernel::stencilLocalSizes1(
    IndexType sizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const common::Grid::BorderType gridBorders[],
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
            IndexType pos[] = { i };

            bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[p], gridSizes, gridBorders, 1 );

            if ( valid )
            {
                cnt++;
            }
        }

        // ToDo: a reflecting boundary might cause two stencil points to be the same

        sizes[ gridPos ] = cnt;
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilLocalSizes2(
    IndexType sizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const common::Grid::BorderType gridBorders[],
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
                IndexType pos[] = { i, j };

                bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[2 * p], gridSizes, gridBorders, 2 );

                if ( valid )
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
    const common::Grid::BorderType gridBorders[],
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
                    IndexType pos[] = { i, j, k };

                    bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[3 * p], gridSizes, gridBorders, 3 );

                    if ( valid )
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
    const common::Grid::BorderType gridBorders[],
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
                        IndexType pos[] = { i, j, k, m};

                        bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[4 * p], gridSizes, gridBorders, 4 );

                        if ( valid )
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    SCAI_REGION( "OpenMP.Stencil.LocalSizes" )


    switch ( nDims ) 
    {
        case 1 : stencilLocalSizes1( sizes, gridSizes, gridDistances, gridBorders, nPoints, stencilNodes );
                 break;

        case 2 : stencilLocalSizes2( sizes, gridSizes, gridDistances, gridBorders, nPoints, stencilNodes );
                 break;

        case 3 : stencilLocalSizes3( sizes, gridSizes, gridDistances, gridBorders, nPoints, stencilNodes );
                 break;

        case 4 : stencilLocalSizes4( sizes, gridSizes, gridDistances, gridBorders, nPoints, stencilNodes );
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int /* stencilOffset */ []  )
{
    #pragma omp parallel for

    for ( IndexType i = 0; i < gridSizes[0]; ++i )
    {
        IndexType gridPos = i * gridDistances[0];

        IndexType ia = csrIA[ gridPos ];

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType pos[] = { i };

            bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[p], gridSizes, gridBorders, 1 );

            if ( !valid )
            {   
                continue;
            }

            csrJA[ ia ] = pos[0] * gridDistances[0];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] /* stencilOffset[] */ )
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
                IndexType pos[] = { i, j };
    
                bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[2 * p], gridSizes, gridBorders, 2 );

                if ( !valid )
                {   
                    continue;
                }

                csrJA[ ia ] = pos[0] * gridDistances[0] + pos[1] * gridDistances[1];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
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
                    IndexType pos[] = { i, j, k };
    
                    bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[3 * p], gridSizes, gridBorders, 3 );

                    if ( !valid )
                    {   
                        continue;
                    }

                    csrJA[ ia ] = pos[0] * gridDistances[0] + pos[1] * gridDistances[1] + pos[2] * gridDistances[2];
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
void OpenMPStencilKernel::stencilLocalCSR4(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
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
                        IndexType pos[] = { i, j, k, m };
    
                        bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[4 * p], gridSizes, gridBorders, 4 );
     
                        if ( !valid )
                        {   
                            continue;
                        }

                        csrJA[ ia ] = pos[0] * gridDistances[0] + pos[1] * gridDistances[1] + pos[2] * gridDistances[2] + pos[3] * gridDistances[3];
                        csrValues[ ia ] = stencilVal[p];

                        ia++;
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.LocalCSR" )

    switch ( nDims ) 
    {
        case 1 : stencilLocalCSR1( csrJA, csrValues, csrIA, gridSizes, gridDistances, gridBorders,
                                   nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilLocalCSR2( csrJA, csrValues, csrIA, gridSizes, gridDistances, gridBorders,
                                   nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilLocalCSR3( csrJA, csrValues, csrIA, gridSizes, gridDistances, gridBorders,
                                   nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilLocalCSR4( csrJA, csrValues, csrIA, gridSizes, gridDistances, gridBorders,
                                   nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilLocalCSR for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

/** Help routine that checks if a certain grid position is local 
 *
 *  @param[in] pos is the position to check
 *  @param[in] lb  are the lower bound in each dimension of the local gri
 *  @param[in] sizes are the local sizes
 *  @param[in] nDims number of dimensions
 */
static inline bool isLocal( IndexType pos[], const IndexType lb[], const IndexType sizes[], const IndexType nDims )
{
    for ( IndexType i = 0; i < nDims; ++i )
    {
        if ( pos[i] < lb[i] )
        {
            return false;
        }
        if ( pos[i] >= lb[i] + sizes[i] )
        {
            return false;
        }
    }
    return true;
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::stencilHaloSizes1(
    IndexType sizes[],
    const IndexType localGridSizes[],
    const IndexType localGridDistances[],
    const IndexType localLB[],
    const IndexType globalGridSizes[],
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    const IndexType nDims = 1;

    // traverse over all points of the local grid

    SCAI_LOG_INFO( logger, "stencilHaloSizes1, local grid size = " << localGridSizes[0] 
                           << ", global grid size = " << globalGridSizes[0] )

    #pragma omp parallel for

    for ( IndexType iLocal = 0; iLocal < localGridSizes[0]; ++iLocal )
    {
        IndexType iGlobal = localLB[0] + iLocal;

        IndexType localIndex = iLocal * localGridDistances[0];

        IndexType cnt = 0;   // check for non-local stencil points but valid in global grid

        // check for each stencil point if it is not local

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType globalPos[] = { iGlobal };

            bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );

            if ( !valid )
            {
                continue;     // stencil point is not available at all
            }

            if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
            {
                continue;     // stencil point is local, so not halo
            }

            cnt ++;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    const IndexType nDims = 2;

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
                IndexType globalPos[] = { iGlobal, jGlobal };

                bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );

                if ( !valid )
                {
                    continue;     // stencil point is not available at all
                }

                if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                {
                    continue;     // stencil point is local, so not halo
                }

                cnt ++;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    const IndexType nDims = 3;

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
                    IndexType globalPos[] = { iGlobal, jGlobal, kGlobal };
                
                    bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );
                
                    if ( !valid )
                    {   
                        continue;     // stencil point is not available at all
                    }
                
                    if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                    {   
                        continue;     // stencil point is local, so not halo
                    }
                
                    cnt ++;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    const IndexType nDims = 4;

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
                        IndexType globalPos[] = { iGlobal, jGlobal, kGlobal, mGlobal };
                    
                        bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );
                    
                        if ( !valid )
                        {   
                            continue;     // stencil point is not available at all
                        }
                    
                        if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                        {   
                            continue;     // stencil point is local, so not halo
                        }
                    
                        cnt ++;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[] )
{
    SCAI_REGION( "OpenMP.Stencil.HaloSizes" )

    switch ( nDims ) 
    {
        case 1 : stencilHaloSizes1( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, globalGridBorders, nPoints, stencilNodes );
                 break;

        case 2 : stencilHaloSizes2( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, globalGridBorders, nPoints, stencilNodes );
                 break;

        case 3 : stencilHaloSizes3( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, globalGridBorders, nPoints, stencilNodes );
                 break;

        case 4 : stencilHaloSizes4( sizes, localGridSizes, localGridDistances, 
                                    localLB, globalGridSizes, globalGridBorders, nPoints, stencilNodes );
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    const IndexType nDims = 1;

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

        // check for each stencil point if it is not local

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType globalPos[] = { iGlobal };

            bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );

            if ( !valid )
            {
                continue;     // stencil point is not available at all
            }

            if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
            {
                continue;     // stencil point is local, so not halo
            }

            csrJA[offset] = globalPos[0] * globalGridDistances[0];
            csrValues[offset] = stencilVal[p];
            ++offset;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    const IndexType nDims = 2;

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

            // check for each stencil point if it is not local

            for ( IndexType p = 0; p < nPoints; ++p )
            {
                IndexType globalPos[] = { iGlobal, jGlobal };

                bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );

                if ( !valid )
                {
                    continue;     // stencil point is not available at all
                }
    
                if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                {
                    continue;     // stencil point is local, so not halo
                }

                csrJA[offset] = globalPos[0] * globalGridDistances[0] + globalPos[1] * globalGridDistances[1];
                csrValues[offset] = stencilVal[p];
                ++offset;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    const IndexType nDims = 3;

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

                // check for each stencil point if it is not local

                for ( IndexType p = 0; p < nPoints; ++p )
                {
                    IndexType globalPos[] = { iGlobal, jGlobal, kGlobal };
                
                    bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );
                
                    if ( !valid )
                    {   
                        continue;     // stencil point is not available at all
                    }
                    
                    if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                    {   
                        continue;     // stencil point is local, so not halo
                    }
                    
                    csrJA[offset] = globalPos[0] * globalGridDistances[0] + globalPos[1] * globalGridDistances[1] +
                                    globalPos[2] * globalGridDistances[2];
                    csrValues[offset] = stencilVal[p];
                    ++offset;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    IndexType nDims = 4;

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

                    // check for each stencil point if it is not local

                    for ( IndexType p = 0; p < nPoints; ++p )
                    {
                        IndexType globalPos[] = { iGlobal, jGlobal, kGlobal, mGlobal };
                    
                        bool valid = common::Grid::getOffsetPos( globalPos, &stencilNodes[nDims * p], globalGridSizes, globalGridBorders, nDims );
                    
                        if ( !valid )
                        {   
                            continue;     // stencil point is not available at all
                        }
                    
                        if ( isLocal( globalPos, localLB, localGridSizes, nDims ) )
                        {   
                            continue;     // stencil point is local, so not halo
                        }
                    
                        csrJA[offset] = globalPos[0] * globalGridDistances[0] + globalPos[1] * globalGridDistances[1] +
                                        globalPos[2] * globalGridDistances[2] + globalPos[3] * globalGridDistances[3];
                        csrValues[offset] = stencilVal[p];
                        ++offset;
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
    const common::Grid::BorderType globalGridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "OpenMP.Stencil.HaloCSR" )

    switch ( nDims ) 
    {
        case 1 : stencilHaloCSR1( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, globalGridBorders,
                                  nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilHaloCSR2( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, globalGridBorders,
                                  nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilHaloCSR3( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, globalGridBorders,
                                  nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilHaloCSR4( csrJA, csrValues, csrIA, localGridSizes, localGridDistances, localLB,
                                  globalGridSizes, globalGridDistances, globalGridBorders,
                                  nPoints, stencilNodes, stencilVal, stencilOffset );
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV1Border" )

    const IndexType nDims = 1;

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
            IndexType pos[] = { i };

            bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[nDims * p], gridSizes, gridBorders, nDims );

            if ( !valid )
            {
                continue;
            }

            IndexType stencilLinearPos = pos[0] * gridDistances[0];

            v += stencilVal[p] * x[ stencilLinearPos ];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    const IndexType nDims = 2;

    SCAI_REGION( "OpenMP.Stencil.GEMV2Border" )

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];

    SCAI_LOG_INFO( logger,  "stencilGEMV2Border on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 )

    #ifdef _OPENMP
    bool outer_parallel = ( i1 - i0 ) > 7;    // only needed if OpenMP enabled
    #endif 

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
                IndexType pos[] = { i, j };

                bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[nDims * p], gridSizes, gridBorders, nDims );

                if ( !valid )
                {
                    continue;
                }

                IndexType stencilLinearPos = pos[0] * gridDistances[0] + pos[1] * gridDistances[1];

                v += stencilVal[p] * x[ stencilLinearPos ];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    SCAI_REGION( "OpenMP.Stencil.GEMV3Border" )

    const IndexType nDims = 3;

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];

    SCAI_LOG_INFO( logger,  "stencilGEMV3Border on " << i0 << " - " << i1 
                             << " x " << j0 << " - " << j1 
                             << " x " << k0 << " - " << k1  )

    #ifdef _OPENMP
    bool outer_parallel = ( i1 - i0 ) > 7;    // only needed if OpenMP enabled
    #endif 

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
                    IndexType pos[] = { i, j, k };

                    bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[nDims * p], gridSizes, gridBorders, nDims );

                    if ( !valid )
                    {
                        continue;
                    }

                    IndexType stencilLinearPos =   pos[0] * gridDistances[0] 
                                                 + pos[1] * gridDistances[1] + pos[2] * gridDistances[2];

                    v += stencilVal[p] * x[ stencilLinearPos ];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int[] )
{
    const IndexType nDims = 4;

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

    #ifdef _OPENMP
    bool outer_parallel = ( i1 - i0 ) > 7;    // only needed if OpenMP enabled
    #endif 

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
                        IndexType pos[] = { i, j, k, m };

                        bool valid = common::Grid::getOffsetPos( pos, &stencilNodes[nDims * p], gridSizes, gridBorders, nDims );
    
                        if ( !valid )
                        {
                            continue;
                        }
    
                        IndexType stencilLinearPos =   pos[0] * gridDistances[0] + pos[1] * gridDistances[1] 
                                                     + pos[2] * gridDistances[2] + pos[3] * gridDistances[3];
    
                        v += stencilVal[p] * x[ stencilLinearPos ];
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
    const common::Grid::BorderType gridBorders[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    switch ( nDims ) 
    {
        case 1 : stencilGEMV1Border( result, alpha, x, gridSizes, gridBounds, gridDistances, gridBorders,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilGEMV2Border( result, alpha, x, gridSizes, gridBounds, gridDistances, gridBorders,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilGEMV3Border( result, alpha, x, gridSizes, gridBounds, gridDistances, gridBorders,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilGEMV4Border( result, alpha, x, gridSizes, gridBounds, gridDistances, gridBorders,
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
    const common::Grid::BorderType gridBorders[],
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

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, gridBorders, nPoints, stencilNodes, stencilVal, stencilOffset );
 
        return;
    }

    // if left boundary, handle it

    if ( width[2 * currentDim] > 0 )
    {
        gridBounds[2 * currentDim] = 0;
        gridBounds[2 * currentDim + 1] = width[2 * currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, gridBorders, nPoints, stencilNodes, stencilVal, stencilOffset );
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

        stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, width, gridDistances, gridBorders, currentDim + 1, 
                           nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // if right boundary, travere it

    if ( width[2 * currentDim + 1] > 0 )
    {
        gridBounds[2 * currentDim] = gridSizes[currentDim] - width[2 * currentDim + 1];
        gridBounds[2 * currentDim + 1] = gridSizes[currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, gridBorders,
                           nPoints, stencilNodes, stencilVal, stencilOffset );
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
    const common::Grid::BorderType gridBorders[],
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

    stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, width, gridDistances, gridBorders, currentDim,
                       nPoints, stencilNodes, stencilVal, stencilOffset );
}

/* --------------------------------------------------------------------------- */

void OpenMPStencilKernel::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry; 

    common::ContextType ctx = common::ContextType::Host;

    SCAI_LOG_DEBUG( logger, "register StencilKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<StencilKernelTrait::stencilLocalSizes>( stencilLocalSizes, ctx, flag );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPStencilKernel::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::ContextType ctx = common::ContextType::Host;

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
