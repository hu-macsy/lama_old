/**
 * @file StorageMethods.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of static routines for matrix storage
 * @author Thomas Brandes
 * @date 27.04.2011
 * $Id$
 */

// hpp
#include <lama/storage/StorageMethods.hpp>

// others
#include <lama/distribution/Distribution.hpp>
#include <lama/distribution/HaloBuilder.hpp>
#include <lama/distribution/Redistributor.hpp>

#include <lama/openmp/OpenMPCSRUtils.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( _StorageMethods::logger, "Storage.Methods" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::localizeCSR(
    LAMAArray<IndexType>& localIA,
    LAMAArray<IndexType>& localJA,
    LAMAArray<ValueType>& localValues,
    const LAMAArray<IndexType>& globalIA,
    const LAMAArray<IndexType>& globalJA,
    const LAMAArray<ValueType>& globalValues,
    const Distribution& rowDist )
{
    LAMA_REGION( "Storage.localizeCSR" )

    const IndexType globalNumRows = globalIA.size() - 1;

    LAMA_ASSERT_EQUAL_ERROR( globalNumRows, rowDist.getGlobalSize() )

    const IndexType localNumRows = rowDist.getLocalSize();

    LAMA_LOG_INFO( logger,
                   "localzeCSR (#rows = global:" << globalNumRows << ", local:" << localNumRows << ", #values( global ) = " << globalJA.size() << ", rowDist = " << rowDist )

    HostReadAccess<IndexType> rGlobalIA( globalIA );

    // localIA will contain at first sizes, later the offsets

    HostWriteOnlyAccess<IndexType> wLocalIA( localIA, localNumRows + 1 );

    // By using the routine local2global of a distribution we can
    // set the local sizes independently

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < localNumRows; i++ )
    {
        IndexType globalI = rowDist.local2global( i );
        wLocalIA[i] = rGlobalIA[globalI + 1] - rGlobalIA[globalI];
    }

    // build now row offset from running sums, set nnz

    IndexType localNumValues = OpenMPCSRUtils::sizes2offsets( wLocalIA.get(), localNumRows );

    LAMA_LOG_DEBUG( logger, "#values( local ) = " << localNumValues )

    // so we can now allocate s of correct size

    HostWriteOnlyAccess<IndexType> wLocalJA( localJA, localNumValues );
    HostWriteOnlyAccess<ValueType> wLocalValues( localValues, localNumValues );

    HostReadAccess<IndexType> rGlobalJA( globalJA );
    HostReadAccess<ValueType> rGlobalValues( globalValues );

    // Due to the availability of offsets for both storages we
    // can copy the global matrix data to the local one independently

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for ( IndexType i = 0; i < localNumRows; i++ )
    {
        IndexType globalI = rowDist.local2global( i );

        IndexType localOffset = wLocalIA[i];

        for ( IndexType jj = rGlobalIA[globalI]; jj < rGlobalIA[globalI + 1]; ++jj )
        {
            wLocalJA[localOffset] = rGlobalJA[jj];
            wLocalValues[localOffset] = rGlobalValues[jj];
            ++localOffset;
        }

        LAMA_ASSERT_DEBUG( localOffset == wLocalIA[i+1],
                           ", localOffset = " << localOffset << ", expected " << wLocalIA[i+1] )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::replicateCSR(
    LAMAArray<IndexType>& globalIA,
    LAMAArray<IndexType>& globalJA,
    LAMAArray<ValueType>& globalValues,
    const LAMAArray<IndexType>& localIA,
    const LAMAArray<IndexType>& localJA,
    const LAMAArray<ValueType>& localValues,
    const Distribution& rowDist )
{
    LAMA_REGION( "Storage.replicateCSR" )

    const IndexType localNumRows = localIA.size() - 1;

    LAMA_ASSERT_EQUAL_ERROR( localNumRows, rowDist.getLocalSize() )

    // the CSR data is not checked for consistency here

    const IndexType globalNumRows = rowDist.getGlobalSize();

    LAMA_LOG_INFO( logger,
                   "replicateCSR (#rows = local:" << localNumRows << ", global: " << globalNumRows << ", #values(local) = " << localJA.size() << ", rowDist = " << rowDist )

    // gather distributed matrix storage into one global storage replicated on all processors

    const Communicator& comm = rowDist.getCommunicator();

    // In a first step we need sizes of all rows, so build it locally before

    LAMAArray<IndexType> localSizes;

    HostReadAccess<IndexType> rLocalIA( localIA );
    HostWriteOnlyAccess<IndexType> wLocalSizes( localSizes, localNumRows );

    OpenMPCSRUtils::offsets2sizes( wLocalSizes.get(), rLocalIA.get(), localNumRows );

    // global IA will contain sizes at first, but later we need the offsets

    HostWriteOnlyAccess<IndexType> wGlobalIA( globalIA, globalNumRows + 1 );

    rowDist.replicate( wGlobalIA.get(), wLocalSizes.get() );

    LAMA_LOG_DEBUG( logger, comm << ": IA is now replicated, build offset " )

    const IndexType globalNumValues = OpenMPCSRUtils::sizes2offsets( wGlobalIA.get(), globalNumRows );

    LAMA_LOG_DEBUG( logger, comm << ": offsets available, #non-zeros (global) = " << globalNumValues )

    // Distribution supports a replicate routine that works on arrays with jagged rows

    HostReadAccess<IndexType> rLocalJA( localJA );
    HostWriteOnlyAccess<IndexType> wGlobalJA( globalJA, globalNumValues );
    rowDist.replicateRagged( wGlobalJA.get(), rLocalJA.get(), wGlobalIA.get() );
    LAMA_LOG_DEBUG( logger, comm << ": JA is now replicated " )

    HostReadAccess<ValueType> rLocalValues( localValues );
    HostWriteOnlyAccess<ValueType> wGlobalValues( globalValues, globalNumValues );
    rowDist.replicateRagged( wGlobalValues.get(), rLocalValues.get(), wGlobalIA.get() );
    LAMA_LOG_DEBUG( logger, comm << ": Values are now replicated " )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::redistributeCSR(
    LAMAArray<IndexType>& targetIA,
    LAMAArray<IndexType>& targetJA,
    LAMAArray<ValueType>& targetValues,
    const LAMAArray<IndexType>& sourceIA,
    const LAMAArray<IndexType>& sourceJA,
    const LAMAArray<ValueType>& sourceValues,
    const Redistributor& redistributor )
{
    LAMA_REGION( "Storage.redistributeCSR" )

    const IndexType sourceNumRows = redistributor.getSourceLocalSize();
    const IndexType targetNumRows = redistributor.getTargetLocalSize();

    LAMA_ASSERT_EQUAL_ERROR( sourceNumRows, sourceIA.size() - 1 )

    LAMA_LOG_INFO( logger,
                   "redistributeCSR, #source rows = " << sourceNumRows << ", #target rows = " << targetNumRows )

    // Redistribution of row sizes, requires size array

    LAMAArray<IndexType> sourceSizes;

    {
        HostWriteOnlyAccess<IndexType> wSourceSizes( sourceSizes, sourceNumRows );
        HostReadAccess<IndexType> rSourceIA( sourceIA );

        OpenMPCSRUtils::offsets2sizes( wSourceSizes.get(), rSourceIA.get(), sourceNumRows );

        HostWriteOnlyAccess<IndexType> wTargetIA( targetIA, targetNumRows + 1 );
    }

    redistributor.redistribute( targetIA, sourceSizes );

    // redistribution does not build ragged row plan by default

    redistributor.buildRowPlans( targetIA, sourceSizes );

    IndexType targetNumValues;

    {
        HostWriteAccess<IndexType> wTargetSizes( targetIA );
        wTargetSizes.resize( targetNumRows + 1 );

        targetNumValues = OpenMPCSRUtils::sizes2offsets( wTargetSizes.get(), targetNumRows );

        LAMA_LOG_INFO( logger,
                       "redistributeCSR: #source values = " << sourceJA.size() << ", #target values = " << targetNumValues )

        // allocate target array with the correct size

        HostWriteOnlyAccess<IndexType> wTargetJA( targetJA, targetNumValues );
        HostWriteOnlyAccess<ValueType> wTargetValues( targetValues, targetNumValues );
    }

    redistributor.redistributeV( targetJA, targetIA, sourceJA, sourceIA );
    redistributor.redistributeV( targetValues, targetIA, sourceValues, sourceIA );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::exchangeHaloCSR(
    LAMAArray<IndexType>& targetIA,
    LAMAArray<IndexType>& targetJA,
    LAMAArray<ValueType>& targetValues,
    const LAMAArray<IndexType>& sourceIA,
    const LAMAArray<IndexType>& sourceJA,
    const LAMAArray<ValueType>& sourceValues,
    const Halo& halo,
    const Communicator& comm )
{
    LAMA_REGION( "Storage.exchangeHaloCSR" )

    // get the number of rows for the new matrix

    const IndexType numRecvRows = halo.getRequiredPlan().totalQuantity();

    LAMAArray<IndexType> sourceSizes;

    // provides indexes: indirection vector needed for sending

    const LAMAArray<IndexType> providesIndexes = halo.getProvidesIndexes();
    const IndexType numSendRows = providesIndexes.size();

    LAMA_LOG_INFO( logger,
                   "exchange halo matrix, #rows to send = " << numSendRows << ", #rows to recv = " << numRecvRows )

    {

        HostReadAccess<IndexType> ia( sourceIA );
        HostReadAccess<IndexType> indexes( providesIndexes );
        HostWriteOnlyAccess<IndexType> sizes( sourceSizes, numSendRows );

        OpenMPCSRUtils::offsets2sizesGather( sizes.get(), ia.get(), indexes, numSendRows );
    }

    {
        // allocate target IA with the right size

        HostWriteOnlyAccess<IndexType> tmpIA( targetIA, numRecvRows + 1 );
    }

    // send the sizes of the rows I will provide

    LAMA_LOG_INFO( logger, "exchange sizes of rows" )

    comm.exchangeByPlan( targetIA, halo.getRequiredPlan(), sourceSizes, halo.getProvidesPlan() );

    // now we build the variable communication plan

    LAMA_LOG_INFO( logger, "exchanged sizes of rows, build vPlans" )

    HostReadAccess<IndexType> sendSizes( sourceSizes );
    HostReadAccess<IndexType> recvSizes( targetIA );

    CommunicationPlan provideV( halo.getProvidesPlan(), sendSizes.get() );
    CommunicationPlan requiredV( halo.getRequiredPlan(), recvSizes.get() );

    const IndexType sendVSize = provideV.totalQuantity();

    LAMA_LOG_INFO( logger, "built vPlans: #send " << sendVSize << ", #recv = " << requiredV.totalQuantity() )

    recvSizes.release();

    LAMA_LOG_INFO( logger, "released read access to recvSizes" )

    // row sizes of target will now become offsets

    {
        HostWriteAccess<IndexType> ia( targetIA );
        ia.resize( numRecvRows + 1 );
        OpenMPCSRUtils::sizes2offsets( ia.get(), numRecvRows );
    }

    // sendJA, sendValues must already be allocated before calling gatherV
    // as its size cannot be determined by arguments

    LAMAArray<IndexType> sendJA( sendVSize );
    LAMAArray<ValueType> sendValues( sendVSize );

    Redistributor::gatherV( sendJA, sourceJA, sourceIA, providesIndexes );
    comm.exchangeByPlan( targetJA, requiredV, sendJA, provideV );

    Redistributor::gatherV( sendValues, sourceValues, sourceIA, providesIndexes );
    comm.exchangeByPlan( targetValues, requiredV, sendValues, provideV );

    // thats it !!
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::splitCSR(
    LAMAArray<IndexType>& localIA,
    LAMAArray<IndexType>& localJA,
    LAMAArray<ValueType>& localValues,
    LAMAArray<IndexType>& haloIA,
    LAMAArray<IndexType>& haloJA,
    LAMAArray<ValueType>& haloValues,
    const LAMAArray<IndexType>& csrIA,
    const LAMAArray<IndexType>& csrJA,
    const LAMAArray<ValueType>& csrValues,
    const Distribution& colDist,
    const Distribution* rowDist )
{
    LAMA_REGION( "Storage.splitCSR" )

    IndexType numRows = csrIA.size() - 1;

    LAMA_LOG_INFO( logger,
                   "splitCSR (#rows = " << numRows << ", #values = " << csrJA.size() << ", colDist = " << colDist )

    if ( rowDist )
    {
        numRows = rowDist->getLocalSize();
    }

    HostReadAccess<IndexType> ia( csrIA );
    HostReadAccess<IndexType> ja( csrJA );

    HostWriteOnlyAccess<IndexType> wLocalIA( localIA, numRows + 1 );
    HostWriteOnlyAccess<IndexType> wHaloIA( haloIA, numRows + 1 );

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType globalI = i;

        if ( rowDist )
        {
            globalI = rowDist->local2global( i );
        }

        IndexType localNonZeros = 0;
        IndexType haloNonZeros = 0;

        for ( IndexType jj = ia[globalI]; jj < ia[globalI + 1]; ++jj )
        {
            if ( colDist.isLocal( ja[jj] ) )
            {
                ++localNonZeros;
            }
            else
            {
                ++haloNonZeros;
            }
        }

        // assignment takes also care for correct first touch

        wLocalIA[i] = localNonZeros;
        wHaloIA[i] = haloNonZeros;
    }

    // build now row offset from running sums, set nnz

    IndexType localNumValues = OpenMPCSRUtils::sizes2offsets( wLocalIA.get(), numRows );
    IndexType haloNumValues = OpenMPCSRUtils::sizes2offsets( wHaloIA.get(), numRows );

    LAMA_LOG_DEBUG( logger,
                    "split: local part has " << localNumValues << " values" << ", halo part has " << haloNumValues << " values" )

    // so we can now allocate s of correct size

    HostWriteOnlyAccess<IndexType> wLocalJA( localJA, localNumValues );
    HostWriteOnlyAccess<ValueType> wLocalValues( localValues, localNumValues );

    HostWriteOnlyAccess<IndexType> wHaloJA( haloJA, haloNumValues );
    HostWriteOnlyAccess<ValueType> wHaloValues( haloValues, haloNumValues );

    HostReadAccess<ValueType> values( csrValues );

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType globalI = i;

        if ( rowDist )
        {
            globalI = rowDist->local2global( i );
        }

        IndexType localOffset = wLocalIA[i];
        IndexType haloOffset = wHaloIA[i];

        LAMA_LOG_TRACE( logger, "fill local row " << i << " from " << localOffset << " to " << wLocalIA[i+1] )

        LAMA_LOG_TRACE( logger, "fill halo row " << i << " from " << haloOffset << " to " << wHaloIA[i+1] )

        for ( IndexType jj = ia[globalI]; jj < ia[globalI + 1]; ++jj )
        {
            const IndexType jLocal = colDist.global2local( ja[jj] );

            if ( jLocal != nIndex )
            {
                // Attention: local gets already local column indexes

                wLocalJA[localOffset] = jLocal;
                wLocalValues[localOffset] = values[jj];

                LAMA_LOG_TRACE( logger,
                                "row " << i << " (global = " << globalI << ": j = " << ja[jj] << " is local " << " local offset = " << localOffset )

                ++localOffset;
            }
            else
            {
                // Attention: halo keeps global column indexes

                wHaloJA[haloOffset] = ja[jj];
                wHaloValues[haloOffset] = values[jj];

                LAMA_LOG_TRACE( logger,
                                "row " << i << " (global = " << globalI << ": j = " << ja[jj] << " is not local " << " halo offset = " << haloOffset )

                ++haloOffset;
            }
        }

        LAMA_ASSERT_DEBUG( localOffset == wLocalIA[i+1],
                           ", localOffset = " << localOffset << ", expected " << wLocalIA[i+1] )
        LAMA_ASSERT_DEBUG( haloOffset == wHaloIA[i+1],
                           ", haloOffset = " << haloOffset << ", expected " << wHaloIA[i+1] )
    }
}

/* -------------------------------------------------------------------------- */

void _StorageMethods::buildHalo(
    class Halo& halo,
    LAMAArray<IndexType>& haloJA,
    IndexType& haloSize,
    const Distribution& colDist )
{
    const IndexType haloNumValues = haloJA.size();

    LAMA_LOG_INFO( logger, "build halo for " << haloNumValues << " non-local indexes" )

    // copy the LAMA array values into a std::vector

    std::vector<IndexType> haloIndexes;
    haloIndexes.reserve( haloNumValues );

    {
        HostReadAccess<IndexType> ja( haloJA );

        for ( IndexType jj = 0; jj < haloNumValues; jj++ )
        {
            haloIndexes.push_back( ja[jj] );
        }
    }

    // sort and then eliminate double elements

    std::sort( haloIndexes.begin(), haloIndexes.end() );

    std::vector<IndexType>::iterator it = std::unique( haloIndexes.begin(), haloIndexes.end() );

    haloIndexes.resize( it - haloIndexes.begin() );

    LAMA_LOG_DEBUG( logger,
                    "Eliminating multiple global indexes, " << haloNumValues << " are shrinked to " << haloIndexes.size() << " values, is halo size" )

    HaloBuilder::build( colDist, haloIndexes, halo );

    LAMA_LOG_DEBUG( logger, "Halo = " << halo )

    haloSize = static_cast<IndexType>( haloIndexes.size() );

    const std::map<IndexType,IndexType>& table = halo.getMap();

    // convert the global indexes to halo indexes

    {
        HostWriteAccess<IndexType> ja( haloJA );

        for ( IndexType jj = 0; jj < haloNumValues; jj++ )
        {
            const IndexType columnIndex = ja[jj];

            const std::map<IndexType,IndexType>::const_iterator elem = table.find( columnIndex );

            LAMA_ASSERT_DEBUG( elem != table.end(), columnIndex << " not found" )

            ja[jj] = elem->second;

            LAMA_ASSERT_DEBUG( ja[jj] < haloSize, "mapped column index out of range " )

            LAMA_LOG_TRACE( logger, "global index " << columnIndex << " is halo index " << ja[jj] )
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageMethods<ValueType>::joinCSR(
    LAMAArray<IndexType>& outIA,
    LAMAArray<IndexType>& outJA,
    LAMAArray<ValueType>& outValues,
    const LAMAArray<IndexType>& localIA,
    const LAMAArray<IndexType>& localJA,
    const LAMAArray<ValueType>& localValues,
    const LAMAArray<IndexType>& haloIA,
    const LAMAArray<IndexType>& haloJA,
    const LAMAArray<ValueType>& haloValues,
    const IndexType numKeepDiagonals )
{
    LAMA_REGION( "Storage.joinCSR" )

    LAMA_ASSERT_EQUAL_ERROR( localIA.size(), haloIA.size() )
    LAMA_ASSERT_EQUAL_ERROR( localJA.size(), localValues.size() )
    LAMA_ASSERT_EQUAL_ERROR( haloJA.size(), haloValues.size() )

    IndexType numRows = localIA.size() - 1;

    LAMA_LOG_INFO( logger,
                   "joinCSRData, #rows = " << numRows << ", local has " << localValues.size() << " elements" << ", halo has " << haloValues.size() << " elements" << ", keep " << numKeepDiagonals << " diagonals " )

    HostWriteOnlyAccess<IndexType> ia( outIA, numRows + 1 );

    HostReadAccess<IndexType> ia1( localIA );
    HostReadAccess<IndexType> ia2( haloIA );

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType nonZeros = 0; // count for each row in parallel

        nonZeros += ia1[i + 1] - ia1[i];
        nonZeros += ia2[i + 1] - ia2[i];

        ia[i] = nonZeros;
    }

    IndexType numValues = OpenMPCSRUtils::sizes2offsets( ia, numRows );

    LAMA_ASSERT_ERROR( numValues == localJA.size() + haloJA.size(), "#non-zero values mismatches" )

    HostWriteOnlyAccess<IndexType> ja( outJA, numValues );
    HostWriteOnlyAccess<ValueType> values( outValues, numValues );

    HostReadAccess<IndexType> ja1( localJA );
    HostReadAccess<ValueType> values1( localValues );

    HostReadAccess<IndexType> ja2( haloJA );
    HostReadAccess<ValueType> values2( haloValues );

    // merging of each row is independent from other rows

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = ia[i];
        IndexType offset1 = ia1[i];
        IndexType offset2 = ia2[i];

        if ( i < numKeepDiagonals )
        {
            if ( offset1 >= ia1[i + 1] )
            {
                LAMA_LOG_FATAL( logger, "no diagonal element for first CSR input data" )
                LAMA_THROWEXCEPTION( "keep diagonal error" )
                // @todo this exception caused segmentation faults when thrown
            }
            ja[offset] = ja1[offset1];
            values[offset++] = values1[offset1++];
        }

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
                LAMA_THROWEXCEPTION( "should not happen here" )
            }
        }

        LAMA_ASSERT_EQUAL_DEBUG( offset1, ia1[i+1] )
        LAMA_ASSERT_EQUAL_DEBUG( offset2, ia2[i+1] )
    }
}

/* -------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT StorageMethods<float> ;
template class LAMA_DLL_IMPORTEXPORT StorageMethods<double> ;

/* -------------------------------------------------------------------------- */

} // namespace LAMA
