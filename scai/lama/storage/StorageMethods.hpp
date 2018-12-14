/**
 * @file StorageMethods.hpp
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

namespace dmemo
{
class Communicator;
class HaloPlan;
class Distribution;
class Redistributor;
}

namespace lama
{

class COMMON_DLL_IMPORTEXPORT _StorageMethods
{
public:

    static void buildHalo(
        class dmemo::HaloPlan& haloPlan,
        hmemo::HArray<IndexType>& haloJA,
        IndexType& haloSize,
        const dmemo::Distribution& colDist );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger ) //!< logger for this matrix format

};

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
        hmemo::HArray<IndexType>& localIA,
        hmemo::HArray<IndexType>& localJA,
        hmemo::HArray<ValueType>& localValues,
        const hmemo::HArray<IndexType>& globalIA,
        const hmemo::HArray<IndexType>& globalJA,
        const hmemo::HArray<ValueType>& globalValues,
        const dmemo::Distribution& rowDist );

    /** Build global CSR storage for row distributed CSR storages.
     *
     *  @param[out]  globalIA, globalJA, globalValues contain the global CSR storage
     *  @param[in]   localIA, localJA, localValues contains the local CSR storage on this proc
     *  @param[in]   rowDist specifies the row distribution of the local CSR stroage
     *
     */

    static void replicateCSR(
        hmemo::HArray<IndexType>& globalIA,
        hmemo::HArray<IndexType>& globalJA,
        hmemo::HArray<ValueType>& globalValues,
        const hmemo::HArray<IndexType>& localIA,
        const hmemo::HArray<IndexType>& localJA,
        const hmemo::HArray<ValueType>& localValues,
        const dmemo::Distribution& rowDist );

    /** Redistribute CSR storages.
     *
     *  @param[out]  targetIA, targetJA, targetValues contain the CSR storage with new distribution
     *  @param[in]   sourceIA, sourceJA, sourceValues contains the CSR storage with original distribution
     *  @param[in]   redistributor specifies the source/target distribution of the local CSR storage
     */

    static void redistributeCSR(
        hmemo::HArray<IndexType>& targetIA,
        hmemo::HArray<IndexType>& targetJA,
        hmemo::HArray<ValueType>& targetValues,
        const hmemo::HArray<IndexType>& sourceIA,
        const hmemo::HArray<IndexType>& sourceJA,
        const hmemo::HArray<ValueType>& sourceValues,
        const dmemo::Redistributor& redistributor );

    /** Exchange rows by halo.
     *
     *  @param[out]  targetIA, targetJA, targetValues contain the new CSR storage
     *  @param[in]   sourceIA, sourceJA, sourceValues contains the original CSR storage
     *  @param[in]   haloPlan used for data exchange
     *  @param[in]   comm     specifies the involved processors
     */

    static void exchangeHaloCSR(
        hmemo::HArray<IndexType>& targetIA,
        hmemo::HArray<IndexType>& targetJA,
        hmemo::HArray<ValueType>& targetValues,
        const hmemo::HArray<IndexType>& sourceIA,
        const hmemo::HArray<IndexType>& sourceJA,
        const hmemo::HArray<ValueType>& sourceValues,
        const dmemo::HaloPlan& haloPlan,
        const dmemo::Communicator& comm );

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
        hmemo::HArray<IndexType>& localIA,
        hmemo::HArray<IndexType>& localJA,
        hmemo::HArray<ValueType>& localValues,
        hmemo::HArray<IndexType>& haloIA,
        hmemo::HArray<IndexType>& haloJA,
        hmemo::HArray<ValueType>& haloValues,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const dmemo::Distribution& colDist,
        const dmemo::Distribution* rowDist );

    /**
     *  Static method that joins rows of two data sets of CSR data.
     *
     *
     *  @param[out] csrIA, csrJA, csrValues will contain the local CSR storage
     *  @param[in]  localIA, localJA, localValues containing the local CSR storage
     *  @param[in]  haloIA, haloJA, haloValues containing the halo CSR storage
     *
     *  The data of one row is sorted according to the column indexes.
     */

    static void joinCSR(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        const hmemo::HArray<IndexType>& localIA,
        const hmemo::HArray<IndexType>& localJA,
        const hmemo::HArray<ValueType>& localValues,
        const hmemo::HArray<IndexType>& haloIA,
        const hmemo::HArray<IndexType>& haloJA,
        const hmemo::HArray<ValueType>& haloValues );
};

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
