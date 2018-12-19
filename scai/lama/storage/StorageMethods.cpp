/**
 * @file StorageMethods.cpp
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
 * @brief Implementation of static routines for matrix storage
 * @author Thomas Brandes
 * @date 27.04.2011
 */

// hpp
#include <scai/lama/storage/StorageMethods.hpp>

// internal scai libraries
#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/HaloExchangePlan.hpp>
#include <scai/dmemo/RedistributePlan.hpp>

#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/utilskernel/TransferUtils.hpp>

#include <scai/common/macros/assert.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <algorithm>

namespace scai
{

using namespace hmemo;
using namespace dmemo;

using sparsekernel::OpenMPCSRUtils;
using sparsekernel::CSRUtils;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( _StorageMethods::logger, "Storage.Methods" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::localizeCSR(
    HArray<IndexType>& localIA,
    HArray<IndexType>& localJA,
    HArray<ValueType>& localValues,
    const HArray<IndexType>& globalIA,
    const HArray<IndexType>& globalJA,
    const HArray<ValueType>& globalValues,
    const Distribution& rowDist )
{
    SCAI_REGION( "Storage.localizeCSR" )
    const IndexType globalNumRows = globalIA.size() - 1;
    SCAI_ASSERT_EQUAL_ERROR( globalNumRows, rowDist.getGlobalSize() )
    const IndexType localNumRows = rowDist.getLocalSize();
    SCAI_LOG_INFO( logger,
                   "localzeCSR (#rows = global:" << globalNumRows << ", local:" << localNumRows << ", #values( global ) = " << globalJA.size() << ", rowDist = " << rowDist )
    ReadAccess<IndexType> rGlobalIA( globalIA );
    // localIA will contain at first sizes, later the offsets
    WriteOnlyAccess<IndexType> wLocalIA( localIA, localNumRows + 1 );
    // By using the routine local2Global of a distribution we can
    // set the local sizes independently
    #pragma omp parallel for 

    for ( IndexType i = 0; i < localNumRows; i++ )
    {
        IndexType globalI = rowDist.local2Global( i );
        wLocalIA[i] = rGlobalIA[globalI + 1] - rGlobalIA[globalI];
    }

    // build now row offset from running sums, set nnz
    IndexType localNumValues = OpenMPCSRUtils::sizes2offsets( wLocalIA.get(), localNumRows );
    SCAI_LOG_DEBUG( logger, "#values( local ) = " << localNumValues )
    // so we can now allocate s of correct size
    WriteOnlyAccess<IndexType> wLocalJA( localJA, localNumValues );
    WriteOnlyAccess<ValueType> wLocalValues( localValues, localNumValues );
    ReadAccess<IndexType> rGlobalJA( globalJA );
    ReadAccess<ValueType> rGlobalValues( globalValues );
    // Due to the availability of offsets for both storages we
    // can copy the global matrix data to the local one independently
    #pragma omp parallel for 

    for ( IndexType i = 0; i < localNumRows; i++ )
    {
        IndexType globalI = rowDist.local2Global( i );
        IndexType localOffset = wLocalIA[i];

        for ( IndexType jj = rGlobalIA[globalI]; jj < rGlobalIA[globalI + 1]; ++jj )
        {
            wLocalJA[localOffset] = rGlobalJA[jj];
            wLocalValues[localOffset] = rGlobalValues[jj];
            ++localOffset;
        }

        SCAI_ASSERT_DEBUG( localOffset == wLocalIA[i + 1],
                           ", localOffset = " << localOffset << ", expected " << wLocalIA[i + 1] )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::replicateCSR(
    HArray<IndexType>& globalIA,
    HArray<IndexType>& globalJA,
    HArray<ValueType>& globalValues,
    const HArray<IndexType>& localIA,
    const HArray<IndexType>& localJA,
    const HArray<ValueType>& localValues,
    const Distribution& rowDist )
{
    SCAI_REGION( "Storage.replicateCSR" )
    const IndexType localNumRows = localIA.size() - 1;
    SCAI_ASSERT_EQUAL_ERROR( localNumRows, rowDist.getLocalSize() )
    // the CSR data is not checked for consistency here
    const IndexType globalNumRows = rowDist.getGlobalSize();
    SCAI_LOG_INFO( logger,
                   "replicateCSR (#rows = local:" << localNumRows << ", global: " << globalNumRows << ", #values(local) = " << localJA.size() << ", rowDist = " << rowDist )
    // gather distributed matrix storage into one global storage replicated on all processors
    const Communicator& comm = rowDist.getCommunicator();
    // In a first step we need sizes of all rows, so build it locally before
    HArray<IndexType> localSizes;
    ReadAccess<IndexType> rLocalIA( localIA );
    WriteOnlyAccess<IndexType> wLocalSizes( localSizes, localNumRows );
    OpenMPCSRUtils::offsets2sizes( wLocalSizes.get(), rLocalIA.get(), localNumRows );
    // global IA will contain sizes at first, but later we need the offsets
    WriteOnlyAccess<IndexType> wGlobalIA( globalIA, globalNumRows + 1 );
    rowDist.replicate( wGlobalIA.get(), wLocalSizes.get() );
    SCAI_LOG_DEBUG( logger, comm << ": IA is now replicated, build offset " )
    const IndexType globalNumValues = OpenMPCSRUtils::sizes2offsets( wGlobalIA.get(), globalNumRows );
    SCAI_LOG_DEBUG( logger, comm << ": offsets available, #non-zeros (global) = " << globalNumValues )
    // Distribution supports a replicate routine that works on arrays with jagged rows
    ReadAccess<IndexType> rLocalJA( localJA );
    WriteOnlyAccess<IndexType> wGlobalJA( globalJA, globalNumValues );
    rowDist.replicateRagged( wGlobalJA.get(), rLocalJA.get(), wGlobalIA.get() );
    SCAI_LOG_DEBUG( logger, comm << ": JA is now replicated " )
    ReadAccess<ValueType> rLocalValues( localValues );
    WriteOnlyAccess<ValueType> wGlobalValues( globalValues, globalNumValues );
    rowDist.replicateRagged( wGlobalValues.get(), rLocalValues.get(), wGlobalIA.get() );
    SCAI_LOG_DEBUG( logger, comm << ": Values are now replicated " )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::redistributeCSR(
    HArray<IndexType>& targetIA,
    HArray<IndexType>& targetJA,
    HArray<ValueType>& targetValues,
    const HArray<IndexType>& sourceIA,
    const HArray<IndexType>& sourceJA,
    const HArray<ValueType>& sourceValues,
    const RedistributePlan& plan )
{
    ContextPtr hostCtx = Context::getHostPtr();

    // Step 1: redistributge sizes (NOT offsets) via redistribution plan

    HArray<IndexType> sourceSizes;
    CSRUtils::offsets2sizes( sourceSizes, sourceIA, hostCtx );
    HArray<IndexType> targetSizes;
    plan.redistribute( targetSizes, sourceSizes );

    // Step 2: build new offset array for target, also needed for exchange of ja and values

    CSRUtils::sizes2offsets( targetIA, targetSizes, hostCtx );
    
    // redistribute ja, and values using variant sizes version

    plan.redistributeV( targetJA, targetIA, sourceJA, sourceIA );
    plan.redistributeV( targetValues, targetIA, sourceValues, sourceIA );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::exchangeHaloCSR(
    HArray<IndexType>& targetIA,
    HArray<IndexType>& targetJA,
    HArray<ValueType>& targetValues,
    const HArray<IndexType>& sourceIA,
    const HArray<IndexType>& sourceJA,
    const HArray<ValueType>& sourceValues,
    const HaloExchangePlan& haloPlan,
    const Communicator& comm )
{
    SCAI_REGION( "Storage.exchangeHaloCSR" )

    ContextPtr loc = Context::getHostPtr();

    HArray<IndexType> sourceSizes; 
    HArray<IndexType> targetSizes;  

    // prepare the sizes of rows to send

    CSRUtils::gatherSizes( sourceSizes, sourceIA, haloPlan.getLocalIndexes(), loc );

    haloPlan.updateHaloDirect( targetSizes, sourceSizes, comm );

    CSRUtils::sizes2offsets( targetIA, targetSizes, loc );

    // now update JA and Values with ragged update

    SCAI_LOG_INFO( logger, "exchanged sizes of rows, build vPlans" )

    // now we build the variable communication plan for JA and Values

    auto provideV = haloPlan.getLocalCommunicationPlan().constructRagged( sourceSizes );
    auto requiredV = haloPlan.getHaloCommunicationPlan().constructRagged( targetSizes );

    const HArray<IndexType>& localIndexes = haloPlan.getLocalIndexes();

    // sendJA, sendValues must already be allocated before calling gatherV
    // as its size cannot be determined by arguments

    const IndexType sendVSize = provideV.totalQuantity();

    HArray<IndexType> sendJA( sendVSize );
    utilskernel::TransferUtils::gatherV( sendJA, sourceJA, sourceIA, localIndexes );
    comm.exchangeByPlan( targetJA, requiredV, sendJA, provideV );

    HArray<ValueType> sendValues( sendVSize );
    utilskernel::TransferUtils::gatherV( sendValues, sourceValues, sourceIA, localIndexes );
    comm.exchangeByPlan( targetValues, requiredV, sendValues, provideV );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::splitCSR(
    HArray<IndexType>& localIA,
    HArray<IndexType>& localJA,
    HArray<ValueType>& localValues,
    HArray<IndexType>& haloIA,
    HArray<IndexType>& haloJA,
    HArray<ValueType>& haloValues,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const Distribution& colDist,
    const Distribution* rowDist )
{
    if ( rowDist )
    {
        if ( rowDist->getGlobalSize() != rowDist->getLocalSize() )
        {
            HArray<IndexType> myIA;
            HArray<IndexType> myJA;
            HArray<ValueType> myValues;
            localizeCSR( myIA, myJA, myValues, csrIA, csrJA, csrValues, *rowDist );
            SCAI_LOG_INFO( logger, "localized before split: ia = " << myIA << ", ja = " << myJA << ", values = " << myValues )
            splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, myIA, myJA, myValues, colDist, NULL );
            return;
        }
    }

    SCAI_REGION( "Storage.splitCSR" )

    ContextPtr ctx = Context::getHostPtr();  // all done on host here

    IndexType numRows = csrIA.size() - 1;

    SCAI_LOG_INFO( logger, "splitCSR( #rows = " << numRows << ", #values = " << csrJA.size()
                            << ", colDist = " << colDist << " ) on " << *ctx )

    HArray<IndexType> isLocalJA;    // contains translation global->local indexes, invalidIndex if not local

    colDist.global2LocalV( isLocalJA, csrJA );

    ReadAccess<IndexType> ia( csrIA, ctx );
    ReadAccess<IndexType> jaLocal( isLocalJA, ctx );
    WriteOnlyAccess<IndexType> wLocalIA( localIA, ctx, numRows + 1 );
    WriteOnlyAccess<IndexType> wHaloIA( haloIA, ctx, numRows + 1 );

    #pragma omp parallel
    {
        SCAI_REGION( "Storage.splitCount" )

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType localNonZeros = 0;
            IndexType haloNonZeros = 0;

            for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
            {
                if ( jaLocal[jj] == invalidIndex )
                {
                    ++haloNonZeros;
                }
                else
                {
                    ++localNonZeros;
                }
            }

            // assignment takes also care for correct first touch
            wLocalIA[i] = localNonZeros;
            wHaloIA[i] = haloNonZeros;
        }
    }

    // build now row offset from running sums, set nnz
    IndexType localNumValues = OpenMPCSRUtils::sizes2offsets( wLocalIA.get(), numRows );
    IndexType haloNumValues = OpenMPCSRUtils::sizes2offsets( wHaloIA.get(), numRows );
    SCAI_LOG_DEBUG( logger,
                    "split: local part has " << localNumValues << " values" << ", halo part has " << haloNumValues << " values" )

    // so we can now allocate s of correct size

    WriteOnlyAccess<IndexType> wLocalJA( localJA, ctx, localNumValues );
    WriteOnlyAccess<ValueType> wLocalValues( localValues, ctx, localNumValues );
    WriteOnlyAccess<IndexType> wHaloJA( haloJA, ctx, haloNumValues );
    WriteOnlyAccess<ValueType> wHaloValues( haloValues, ctx, haloNumValues );

    ReadAccess<IndexType> jaGlobal( csrJA, ctx );
    ReadAccess<ValueType> values( csrValues, ctx );

    #pragma omp parallel
    {
        SCAI_REGION( "Storage.splitTransfer" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType localOffset = wLocalIA[i];
            IndexType haloOffset = wHaloIA[i];

            SCAI_LOG_TRACE( logger, "fill local row " << i << " from " << localOffset << " to " << wLocalIA[i + 1] )
            SCAI_LOG_TRACE( logger, "fill halo row " << i << " from " << haloOffset << " to " << wHaloIA[i + 1] )

            for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
            {
                const IndexType jLocal = jaLocal[jj];

                if ( jLocal != invalidIndex )
                {
                    // local column index, will be used in local JA

                    wLocalJA[localOffset] = jLocal;
                    wLocalValues[localOffset] = values[jj];

                    SCAI_LOG_TRACE( logger,
                                    "row " << i << " : j = " << jaGlobal[jj]
                                    << " is local " << " local offset = " << localOffset )

                    ++localOffset;
                }
                else
                {
                    // non-local index, global index remains in halo JA

                    wHaloJA[haloOffset] = jaGlobal[jj];    // global -> halo translation is done later
                    wHaloValues[haloOffset] = values[jj];

                    SCAI_LOG_TRACE( logger,
                                    "row " << i << " (: j = " << jaGlobal[jj]
                                    << " is not local " << " halo offset = " << haloOffset )
                    ++haloOffset;
                }
            }

            SCAI_ASSERT_EQ_DEBUG( localOffset, wLocalIA[i + 1], "serious mismatch to previsous counting" )
            SCAI_ASSERT_EQ_DEBUG( haloOffset, wHaloIA[i + 1], "serious mitmatch to previous counting" )
        }
    }
}

/* -------------------------------------------------------------------------- */

void _StorageMethods::buildHaloExchangePlan(
    HaloExchangePlan& haloPlan,
    HArray<IndexType>& haloJA,
    const Distribution& colDist )
{
    SCAI_REGION( "Storage.buildHaloExchangePlan" )

    const IndexType haloNumValues = haloJA.size();

    SCAI_LOG_INFO( logger, "build halo for " << haloNumValues << " non-local indexes" )

    bool elimDoubleHere = true;

    if ( elimDoubleHere )
    {
        // eliminate double elements in haloJA is more efficient here:
        //  - sort + unique faster than using a map 

        HArray<IndexType> requiredIndexes( haloJA );

        {
            SCAI_REGION( "Storage.elimDoubleRequired" )
            auto wIndexes = hostWriteAccess( requiredIndexes );
            std::sort( wIndexes.begin(), wIndexes.end() );
            auto it = std::unique( wIndexes.begin(), wIndexes.end() );
            wIndexes.resize( it - wIndexes.begin() );
        }

        haloPlan = haloExchangePlanByRequiredIndexes( requiredIndexes, colDist, false );
    }
    else
    {
        haloPlan = haloExchangePlanByRequiredIndexes( haloJA, colDist, true );
    }

    SCAI_LOG_DEBUG( logger, "Halo plan = " << haloPlan )

    // translate the non-local - global indexes to halo indexes, in-place

    haloPlan.global2HaloV( haloJA, haloJA );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::joinCSR(
    HArray<IndexType>& outIA,
    HArray<IndexType>& outJA,
    HArray<ValueType>& outValues,
    const HArray<IndexType>& localIA,
    const HArray<IndexType>& localJA,
    const HArray<ValueType>& localValues,
    const HArray<IndexType>& haloIA,
    const HArray<IndexType>& haloJA,
    const HArray<ValueType>& haloValues )
{
    SCAI_REGION( "Storage.joinCSR" )
    SCAI_ASSERT_EQUAL_ERROR( localIA.size(), haloIA.size() )
    SCAI_ASSERT_EQUAL_ERROR( localJA.size(), localValues.size() )
    SCAI_ASSERT_EQUAL_ERROR( haloJA.size(), haloValues.size() )
    IndexType numRows = localIA.size() - 1;
    SCAI_LOG_INFO( logger,
                   "joinCSRData, #rows = " << numRows << ", local has " << localValues.size() << " elements" 
                   << ", halo has " << haloValues.size() << " elements" )
    WriteOnlyAccess<IndexType> ia( outIA, numRows + 1 );
    ReadAccess<IndexType> ia1( localIA );
    ReadAccess<IndexType> ia2( haloIA );
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType nonZeros = 0; // count for each row in parallel
        nonZeros += ia1[i + 1] - ia1[i];
        nonZeros += ia2[i + 1] - ia2[i];
        ia[i] = nonZeros;
    }

    IndexType numValues = OpenMPCSRUtils::sizes2offsets( ia, numRows );
    SCAI_ASSERT_ERROR( numValues == localJA.size() + haloJA.size(), "#non-zero values mismatches" )
    WriteOnlyAccess<IndexType> ja( outJA, numValues );
    WriteOnlyAccess<ValueType> values( outValues, numValues );
    ReadAccess<IndexType> ja1( localJA );
    ReadAccess<ValueType> values1( localValues );
    ReadAccess<IndexType> ja2( haloJA );
    ReadAccess<ValueType> values2( haloValues );
    // merging of each row is independent from other rows
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = ia[i];
        IndexType offset1 = ia1[i];
        IndexType offset2 = ia2[i];

        // merge other data by order of column indexes

        for ( IndexType k = offset; k < ia[i + 1]; ++k )
        {
            // merge the row data of local and halo
            bool isNext1 = offset1 < ia1[i + 1];
            bool isNext2 = offset2 < ia2[i + 1]; // if same type

            if ( isNext1 && isNext2 )
            {
                // both arrays not empty, take the one with smaller global index
                isNext1 = ja1[offset1] < ja2[offset2];
                isNext2 = !isNext1;
            }

            if ( isNext1 )
            {
                ja[k] = ja1[offset1];
                values[k] = values1[offset1++];
            }
            else if ( isNext2 )
            {
                ja[k] = ja2[offset2];
                values[k] = values2[offset2++];
            }
            else
            {
                COMMON_THROWEXCEPTION( "should not happen here" )
            }
        }

        SCAI_ASSERT_EQUAL_DEBUG( offset1, ia1[i + 1] )
        SCAI_ASSERT_EQUAL_DEBUG( offset2, ia2[i + 1] )
    }
}

/* -------------------------------------------------------------------------- */

SCAI_COMMON_INST_CLASS( StorageMethods, SCAI_NUMERIC_TYPES_HOST )

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
