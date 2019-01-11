/**
 * @file TransferUtils.hpp
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
 * @brief Transer of elements in heterogeneuous arrays via indirect addressing (gather/scatter)
 * @author Thomas Brandes
 * @date 08.10.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library

// internal scai libraries

#include <scai/hmemo/HArray.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tracing.hpp>

#include <memory>

namespace scai
{

namespace utilskernel
{

/** Static class of transfer routines of HArray data
 *
 */
class COMMON_DLL_IMPORTEXPORT TransferUtils: public scai::common::Printable
{
public:

    /** 
     *  @brief Gather of rows in a dense matrix, each row has n entries
     *
     *  @param[out] targetArray will contain the gathered rows, n * sourceIndexes.size()
     *  @param[in]  sourceArray contains all rows
     *  @param[in]  sourceIndexes specfies the selected rows 
     *  @param[in]  n is number of entries in one row
     */
    template<typename ValueType>
    static void gatherN(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceIndexes,
        const IndexType n );

    template<typename ValueType>
    static void scatterN(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetIndexes,
        const hmemo::HArray<ValueType>& sourceArray,
        const IndexType n );

    /**
     *  Gather but an individual size for each row
     */
    template<typename ValueType>
    static void gatherV(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceOffsets,
        const hmemo::HArray<IndexType>& sourceIndexes );

    template<typename ValueType>
    static void scatterV(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetOffsets,
        const hmemo::HArray<IndexType>& targetIndexes,
        const hmemo::HArray<ValueType>& sourceArray );

    template<typename ValueType>
    static void copy(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetIndexes,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceIndexes );

    template<typename ValueType>
    static void copyN(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetIndexes,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceIndexes,
        IndexType n );

    template<typename ValueType>
    static void copyV(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetOffsets,
        const hmemo::HArray<IndexType>& targetIndexes,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceOffsets,
        const hmemo::HArray<IndexType>& sourceIndexes );
};

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::gatherN(

    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceIndexes,
    const IndexType n )
{
    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    hmemo::WriteAccess<ValueType> target( targetArray, loc );
    hmemo::ReadAccess<ValueType> source( sourceArray, loc );
    hmemo::ReadAccess<IndexType> indexes( sourceIndexes, loc );
    #pragma omp parallel for

    for ( IndexType i = 0; i < indexes.size(); i++ )
    {
//        SCAI_LOG_DEBUG( logger,
//                        "targetN[" << i << "] = sourceN[" << indexes[i] << "] = " << source[indexes[i] * n] << " ..." )

        for ( IndexType j = 0; j < n; j++ )
        {
            target[i * n + j] = source[indexes[i] * n + j];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::gatherV(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceOffsets,
    const hmemo::HArray<IndexType>& sourceIndexes )
{
    SCAI_REGION( "Transfer.gatherV" )

    const IndexType n = sourceIndexes.size();
    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    hmemo::WriteAccess<ValueType> wTargetArray( targetArray, loc );
    hmemo::ReadAccess<ValueType> rSourceArray( sourceArray, loc );
    hmemo::ReadAccess<IndexType> rSourceOffsets( sourceOffsets, loc );
    hmemo::ReadAccess<IndexType> rSourceIndexes( sourceIndexes, loc );
    // Note: we have no target offsets array

    IndexType targetSize = targetArray.size();
    IndexType targetOffset = 0;

    for ( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = rSourceIndexes[ii];

        SCAI_ASSERT_LE_DEBUG( targetOffset + ( rSourceOffsets[i + 1] - rSourceOffsets[i] ), targetSize,
                              "target array too small" )

        for ( IndexType j = rSourceOffsets[i]; j < rSourceOffsets[i + 1]; ++j )
        {
            wTargetArray[targetOffset++] = rSourceArray[j];
        }
    }

    SCAI_ASSERT_LE_ERROR( targetOffset, targetSize, "target array was too small" )

    targetArray.resize( targetOffset );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::scatterN(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetIndexes,
    const hmemo::HArray<ValueType>& sourceArray,
    const IndexType n )
{
    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    hmemo::WriteAccess<ValueType> target( targetArray, loc );
    hmemo::ReadAccess<IndexType> indexes( targetIndexes, loc );
    hmemo::ReadAccess<ValueType> source( sourceArray, loc );
    #pragma omp parallel for

    for ( IndexType i = 0; i < indexes.size(); i++ )
    {
//        SCAI_LOG_DEBUG( logger,
//                        "targetN[" << indexes[i] << "] = sourceN[" << i << "] = " << source[i * n] << " ..." )

        for ( IndexType j = 0; j < n; j++ )
        {
            target[indexes[i] * n + j] = source[i * n + j];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::scatterV(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetOffsets,
    const hmemo::HArray<IndexType>& targetIndexes,
    const hmemo::HArray<ValueType>& sourceArray )
{
    SCAI_REGION( "Transfer.scatterV" )

    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    const IndexType n = targetIndexes.size();
    hmemo::WriteAccess<ValueType> wTargetArray( targetArray, loc );
    hmemo::ReadAccess<IndexType> rTargetOffsets( targetOffsets, loc );
    hmemo::ReadAccess<IndexType> rTargetIndexes( targetIndexes, loc );
    hmemo::ReadAccess<ValueType> rSourceArray( sourceArray, loc );
    // Note: we have no source offsets array, no parallelization possible
    IndexType sourceOffset = 0;

    for ( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = rTargetIndexes[ii];

        for ( IndexType j = rTargetOffsets[i]; j < rTargetOffsets[i + 1]; ++j )
        {
            wTargetArray[j] = rSourceArray[sourceOffset++];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::copy(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetIndexes,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceIndexes )
{
    SCAI_REGION( "Transfer.copy" )

    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    hmemo::WriteAccess<ValueType> target( targetArray, loc );
    hmemo::ReadAccess<ValueType> source( sourceArray, loc );
    hmemo::ReadAccess<IndexType> tindexes( targetIndexes, loc );
    hmemo::ReadAccess<IndexType> sindexes( sourceIndexes, loc );
    SCAI_ASSERT_ERROR( tindexes.size() == sindexes.size(), "index size mismatch" )

    for ( IndexType i = 0; i < tindexes.size(); i++ )
    {
//        SCAI_LOG_DEBUG( logger,
//                        "target[" << tindexes[i] << "] = source[" << sindexes[i] << "] = " << source[ sindexes[i] ] )
        target[tindexes[i]] = source[sindexes[i]];
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::copyN(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetIndexes,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceIndexes,
    IndexType n )
{
    SCAI_REGION( "Transfer.copyN" )

    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    hmemo::WriteAccess<ValueType> target( targetArray, loc );
    hmemo::ReadAccess<ValueType> source( sourceArray, loc );
    hmemo::ReadAccess<IndexType> tindexes( targetIndexes, loc );
    hmemo::ReadAccess<IndexType> sindexes( sourceIndexes, loc );
    SCAI_ASSERT_ERROR( tindexes.size() == sindexes.size(), "index size mismatch" )
    #pragma omp parallel for

    for ( IndexType i = 0; i < tindexes.size(); i++ )
    {
//        SCAI_LOG_DEBUG( logger,
//                        "targetN[" << tindexes[i] << "] = sourceN[" << sindexes[i] << "] = " << source[ sindexes[i] * n ] << " ..." )

        for ( IndexType j = 0; j < n; j++ )
        {
            target[tindexes[i] * n + j] = source[sindexes[i] * n + j];
        }
    }
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void TransferUtils::copyV(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetOffsets,
    const hmemo::HArray<IndexType>& targetIndexes,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceOffsets,
    const hmemo::HArray<IndexType>& sourceIndexes )
{
    SCAI_REGION( "Transfer.copyV" )

    SCAI_ASSERT_EQ_ERROR( targetIndexes.size(), sourceIndexes.size(), "size mismatch" )
    hmemo::ContextPtr loc = hmemo::Context::getHostPtr();
    const IndexType n = targetIndexes.size();
    hmemo::WriteAccess<ValueType> wTargetArray( targetArray, loc );
    hmemo::ReadAccess<IndexType> rTargetOffsets( targetOffsets, loc );
    hmemo::ReadAccess<IndexType> rTargetIndexes( targetIndexes, loc );
    hmemo::ReadAccess<ValueType> rSourceArray( sourceArray, loc );
    hmemo::ReadAccess<IndexType> rSourceOffsets( sourceOffsets, loc );
    hmemo::ReadAccess<IndexType> rSourceIndexes( sourceIndexes, loc );

    for ( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType sourceI = rSourceIndexes[ii];
        IndexType targetI = rTargetIndexes[ii];
        IndexType k = rTargetOffsets[targetI];

        for ( IndexType j = rSourceOffsets[sourceI]; j < rSourceOffsets[sourceI + 1]; ++j )
        {
            wTargetArray[k] = rSourceArray[j];
            ++k;
        }

        SCAI_ASSERT_EQ_DEBUG( k, rTargetOffsets[targetI + 1], "size mismatch" )
    }
}

} /* end namespace utilskernel */

} /* end namespace scai */
