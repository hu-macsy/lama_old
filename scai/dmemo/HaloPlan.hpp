/**
 * @file HaloPlan.hpp
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
 * @brief HaloPlan.hpp
 * @author Thomas Brandes
 * @date 23.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/dmemo/CommunicationPlan.hpp>
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/BinaryOp.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <map>

namespace scai
{

namespace dmemo
{

/** 
 *  A HaloPlan is an object that describes the communication pattern for the data exchange of
 *  local entries with incarnation in the halo part on other processors that require the actual
 *  value.
 *
 *  The halo on each processor itself is an array  from 0 to haloSize where each of its entries
 *  stands for a non-local entry of another processor where the entry is specified by its global
 *  index required by this processor. The halo plan itself contains the communication plans for updating 
 *  the halo with actual values from the owning processors that provide the data.
 *
 */

class COMMON_DLL_IMPORTEXPORT HaloPlan: public common::Printable
{
public:

    /** Constructor of a new 'empty' halo. */

    HaloPlan();

    /** Constructor 
     *
     *  @param[in] requiredIndexes array of global indexes, already sorted by processors
     *  @param[in] providesIndexes array of local indexes, sorted as required by other processors
     *  @param[in] requiredPlan  communication plan for receiving data according to the required indexes
     *  @param[in] providesPlan  communication plan for sending data according to the providesIndexes
     *
     *  - requiredIndexes.size() == requiredPlan.totalQuantitity()
     *  - providesIndexes.size() == providesPlan.totalQuantitity()
     *  - requiredPlan == transpose( providesPlan )
     *
     *  Note: for halo updates itself the required indexes are not needed but they specify the 
     *        relation between the halo indexes and the global indexes.
     */
    HaloPlan( hmemo::HArray<IndexType> requiredIndexes,
              hmemo::HArray<IndexType> providesIndexes,
              CommunicationPlan requiredPlan, 
              CommunicationPlan providesPlan );

    /** Split up */

    void splitUp( hmemo::HArray<IndexType>& requiredIndexes,
                  hmemo::HArray<IndexType>& providesIndexes,
                  CommunicationPlan& requiredPlan,
                  CommunicationPlan& providesPlan );

    /** Copy constructor. */

    HaloPlan( const HaloPlan& halo ) = default;
    HaloPlan( HaloPlan&& halo ) = default;

    virtual ~HaloPlan();

    /** Clear the halo plan, but do not free the allocated data */

    void clear();

    /** Purge, same as clear but also frees the allocated data */

    void purge();

    /** Swap member variables with another halo plan. */

    void swap( HaloPlan& other );

    /** 
     *  This method updates the halo (non-local values) with the required values from the source array.
     *
     *  @param[out] haloArray will contain the actual values of the source array for the required indexes
     *  @param[in]  sourceArray  provides the values for update of the halo
     *  @param[in]  comm         is used for the exchange of the data          
     *
     *  Note: each processor that provides or requires data must call this routine.
     */
    template<typename ValueType>
    void updateHalo( 
        hmemo::HArray<ValueType>& haloArray, 
        const hmemo::HArray<ValueType>& sourceArray, 
        const Communicator& comm ) const;

    /**
     *  same as updateHalo but provides a temporary array that takes the send buffer
     *
     *  Optimized version that avoids the allocation of an additional array for the send buffer.
     */
    template<typename ValueType>
    void updateHalo( 
        hmemo::HArray<ValueType>& haloArray, 
        const hmemo::HArray<ValueType>& sourceArray, 
        const Communicator& comm,
        hmemo::HArray<ValueType>& tmp ) const;

    /** 
     *  Update halo data directly with the send data that is the gathered local data.
     *
     *  Splitting update in gather and send is required when communication is overlapped
     *  with an update of the local data, so gathering of data must be done before.
     */
    template<typename ValueType>
    void updateHaloDirect( 
        hmemo::HArray<ValueType>& haloArray, 
        const hmemo::HArray<ValueType>& sendArray,
        const Communicator& comm ) const;

    /**
     *  Asynchronous version
     */
    template<typename ValueType>
    tasking::SyncToken* updateHaloAsync( 
        hmemo::HArray<ValueType>& haloArray, 
        const hmemo::HArray<ValueType>& sourceArray, 
        const Communicator& comm ) const;
    /** 
     *  This method takes the halo (non-local values) to update the source array.
     */
    template<typename ValueType>
    void updateByHalo( hmemo::HArray<ValueType>& sourceArray, 
                       const hmemo::HArray<ValueType>& haloArray, 
                       common::BinaryOp op, const Communicator& comm ) const;

    HaloPlan& operator=( const HaloPlan& other ) = default;
    HaloPlan& operator=( HaloPlan&& other ) = default;

    inline const CommunicationPlan& getHaloCommunicationPlan() const;

    inline const CommunicationPlan& getLocalCommunicationPlan() const;

    inline IndexType global2halo( const IndexType globalIndex ) const;

    inline const hmemo::HArray<IndexType>& getProvidesIndexes() const;

    /**
     * Getter for the array that contains the translation from halo indexes to global indexes.
     */
    inline const hmemo::HArray<IndexType>& getHalo2GlobalIndexes() const;

    /** Query the size for a halo to be allocated */

    inline IndexType getHaloSize() const;

    /** If a halo is empty, no communication is needed for this partition.
     Be careful: getHaloSize() == 0 implies that no indexes are required
     but it might be possible that this partition has to provide values
     */

    inline bool isEmpty() const;

    virtual void writeAt( std::ostream& stream ) const;

    /** 
     *  @brief This method translates a set of halo indexes to the corresponding global indexes.
     *
     *  @param[out] globalIndexes is the array with the indexes that are translated
     *  @param[in]  haloIndexes is the array with halo indexes in range 0, .., haloSize() - 1
     *
     *  Note: alias of localIndexes and haloIndexes is supported
     */
    void halo2GlobalV( hmemo::HArray<IndexType>& globalIndexes, const hmemo::HArray<IndexType>& haloIndexes ) const;

    /** This method translates global indexes into local indexes 
     *
     *  @param[in,out] indexes is the array with global indexes that are translated halo indexes
     *
     *  The array must only contain required indexes used to build the halo.
     */
    void global2HaloV( hmemo::HArray<IndexType>& haloIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const;

    /** Constructor function as static version for more convenient logging */

    static HaloPlan constructByRequiredIndexes( const hmemo::HArray<IndexType>& requiredIndexes, const class Distribution& distribution );

private:

    // Indexes for required values and values to provide are stored in HArrays
    // so they might be used in different contexts, especially also on GPU

    hmemo::HArray<IndexType> mRequiredIndexes;   // non-local global indexes required
    hmemo::HArray<IndexType> mProvidesIndexes;   // local indexes for provided data

    CommunicationPlan mRequiredPlan;
    CommunicationPlan mProvidesPlan;   // is the transpose of mRequiredPlan

    // This is an additional data structure to map the required indexes to halo indexes
    // mGlobal2Halo[mRequiredIndexes[i]] = i

    std::map<IndexType, IndexType> mGlobal2Halo;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

const CommunicationPlan& HaloPlan::getHaloCommunicationPlan() const
{
    return mRequiredPlan;
}

const CommunicationPlan& HaloPlan::getLocalCommunicationPlan() const
{
    return mProvidesPlan;
}

const hmemo::HArray<IndexType>& HaloPlan::getProvidesIndexes() const
{
    return mProvidesIndexes;
}

const hmemo::HArray<IndexType>& HaloPlan::getHalo2GlobalIndexes() const
{
    return mRequiredIndexes;
}

IndexType HaloPlan::global2halo( const IndexType globalIndex ) const
{
    const std::map<IndexType, IndexType>::const_iterator elem = mGlobal2Halo.find( globalIndex );

    if ( elem == mGlobal2Halo.end() )
    {
        return invalidIndex;
    }

    return elem->second;
}

IndexType HaloPlan::getHaloSize() const
{
    // Note: is also same as mRequiredIndexes.size()

    return mRequiredPlan.totalQuantity();
}

bool HaloPlan::isEmpty() const
{
    return ( mRequiredPlan.totalQuantity() == 0 ) && ( mProvidesPlan.totalQuantity() == 0 );
}

/** Build a halo by an array of (unsorted) required indexes 
 *
 *  @param[in]  requiredIndexes are global indexes for required values from other processors
 *  @param[in]  distribution is the mapping to find owners and local indexes 
 *  @returns    the HaloPlan object that contains (sorted) required and provides indexes and exchange plans
 *
 *  Note: requiredIndexes should not contain global indexes that are owned by this processor (isLocal)
 *        and there should be no double values in it to reduce communication volume.
 */
inline HaloPlan haloPlanByRequiredIndexes( const hmemo::HArray<IndexType>& requiredIndexes, const class Distribution& distribution )
{
    return HaloPlan::constructByRequiredIndexes( requiredIndexes, distribution );
}

} /* end namespace dmemo */

} /* end namespace scai */
