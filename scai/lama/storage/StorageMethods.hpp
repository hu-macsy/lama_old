/**
 * @file StorageMethods.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Class providing static routines for matrix storage
 * @author Thomas Brandes
 * @date 27.04.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

using hmemo::HArray;

namespace lama
{

class COMMON_DLL_IMPORTEXPORT _StorageMethods
{
public:

    static void buildHalo(
        class Halo& halo,
        HArray<IndexType>& haloJA,
        IndexType& haloSize,
        const class Distribution& colDist );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format

}    ;

/* -------------------------------------------------------------------------- */

/** This class provides static utility methods for splitting matrix storage into a
 *  local and a halo part.
 *
 *  Due to a column distribution the storage is divided into a local part (having
 *  the local columns) and a halo part (for the non-local columns). Furthermore,
 *  it builds the halo for exchanging the non-local values between processors.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT StorageMethods: public _StorageMethods
{
public:

    /** Localize CSR storage for row distribution.
     *
     *  @param[out] localIA, localJA, localValues will contain the local CSR storage
     *  @param[in]  globalIA, globalJA, globalValues contain the global CSR storage
     *  @param[in]  rowDist distribution so only local rows of CSR storage are taken
     */

    static void localizeCSR(
        HArray<IndexType>& localIA,
        HArray<IndexType>& localJA,
        HArray<ValueType>& localValues,
        const HArray<IndexType>& globalIA,
        const HArray<IndexType>& globalJA,
        const HArray<ValueType>& globalValues,
        const class Distribution& rowDist );

    /** Build global CSR storage for row distributed CSR storages.
     *
     *  @param[out]  globalIA, globalJA, globalValues contain the global CSR storage
     *  @param[in]   localIA, localJA, localValues contains the local CSR storage on this proc
     *  @param[in]   rowDist specifies the row distribution of the local CSR stroage
     *
     */

    static void replicateCSR(
        HArray<IndexType>& globalIA,
        HArray<IndexType>& globalJA,
        HArray<ValueType>& globalValues,
        const HArray<IndexType>& localIA,
        const HArray<IndexType>& localJA,
        const HArray<ValueType>& localValues,
        const class Distribution& rowDist );

    /** Redistribute CSR storages.
     *
     *  @param[out]  targetIA, targetJA, targetValues contain the CSR storage with new distribution
     *  @param[in]   sourceIA, sourceJA, sourceValues contains the CSR storage with original distribution
     *  @param[in]   redistributor specifies the source/target distribution of the local CSR storage
     */

    static void redistributeCSR(
        HArray<IndexType>& targetIA,
        HArray<IndexType>& targetJA,
        HArray<ValueType>& targetValues,
        const HArray<IndexType>& sourceIA,
        const HArray<IndexType>& sourceJA,
        const HArray<ValueType>& sourceValues,
        const class Redistributor& redistributor );

    /** Exchange rows by halo.
     *
     *  TODO[doxy] Complete Description.
     *
     *  @param[out]  targetIA, targetJA, targetValues contain the new CSR storage
     *  @param[in]   sourceIA, sourceJA, sourceValues contains the original CSR storage
     *  @param[in]   halo
     *  @param[in]   comm
     */

    static void exchangeHaloCSR(
        HArray<IndexType>& targetIA,
        HArray<IndexType>& targetJA,
        HArray<ValueType>& targetValues,
        const HArray<IndexType>& sourceIA,
        const HArray<IndexType>& sourceJA,
        const HArray<ValueType>& sourceValues,
        const class Halo& halo,
        const class Communicator& comm );

    /** Splitting CSR storage.
     *
     *  @param[out] localIA, localJA, localValues will contain the local CSR storage
     *  @param[out] haloIA, haloJA, haloValues will contain the halo CSR storage
     *  @param[in]   csrIA, csrJA, csrValues is the storage to split
     *  @param[in]   colDist column distribution used for splitting
     *  @param[in]   rowDist optional row distribution so only local rows of CSR storage are taken
     *
     *  The halo CSR storage will still contain the global indexes.
     */

    static void splitCSR(
        HArray<IndexType>& localIA,
        HArray<IndexType>& localJA,
        HArray<ValueType>& localValues,
        HArray<IndexType>& haloIA,
        HArray<IndexType>& haloJA,
        HArray<ValueType>& haloValues,
        const HArray<IndexType>& csrIA,
        const HArray<IndexType>& csrJA,
        const HArray<ValueType>& csrValues,
        const class Distribution& colDist,
        const class Distribution* rowDist );

    /**
     *  Static method that joins rows of two data sets of CSR data.
     *
     *
     *  @param[out] csrIA, csrJA, csrValues will contain the local CSR storage
     *  @param[in]  localIA, localJA, localValues containing the local CSR storage
     *  @param[in]  haloIA, haloJA, haloValues containing the halo CSR storage
     *  @param[in]  numKeepDiagonals number of rows where diagonal element of local should be taken at first
     *
     *  The data of one row is sorted according to the column indexes.
     */

    static void joinCSR(
        HArray<IndexType>& csrIA,
        HArray<IndexType>& csrJA,
        HArray<ValueType>& csrValues,
        const HArray<IndexType>& localIA,
        const HArray<IndexType>& localJA,
        const HArray<ValueType>& localValues,
        const HArray<IndexType>& haloIA,
        const HArray<IndexType>& haloJA,
        const HArray<ValueType>& haloValues,
        const IndexType numKeepDiagonals );
};

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
