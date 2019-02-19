/**
 * @file HaloExchangePlan.hpp
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
 * @brief HaloExchangePlan.hpp
 * @author Thomas Brandes
 * @date 17.12.2018
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
 *  A HaloExchangePlan is an object that describes the communication pattern for the data exchange of
 *  local entries with incarnation in the halo part on other processors that require the actual
 *  value.
 *
 *  The halo on each processor itself is an array  from 0 to haloSize where each of its entries
 *  stands for a non-local entry of another processor where the entry is specified by its global
 *  index required by this processor. The halo plan itself contains the communication plans for updating 
 *  the halo with actual values from the owning processors that provide the data.
 *
 */
class COMMON_DLL_IMPORTEXPORT HaloExchangePlan: public common::Printable
{
public:

    /** Constructor of a new 'empty' halo. */

    HaloExchangePlan();

    /** Constructor 
     *
     *  @param[in] halo2GlobalIndexes array of global indexes for which halo entries stand
     *  @param[in] localIndexes array of local indexes that have a correspondant halo entry on other processors
     *  @param[in] haloCommPlan  communication plan for exchange of data in the halo array
     *  @param[in] localCommPlan communication plan for exchange of data in the local array
     *  @param[in] global2HaloMap map of global indexes to halo indexes
     *
     *  - halo2GlobalIndexes.size() == haloCommPlan.totalQuantitity()
     *  - localIndexes.size() == localCommPlan.totalQuantitity()
     *  - haloCommPlan == transpose( localCommPlan )
     *
     *  Notes: 
     *
     *  - halo2GlobalIndexes must be sorted by the owners of the global indexes
     *  - the local indexes are sorted so that communication for each processor is contiguous
     *  - for halo updates itself the halo-global indexes are not needed.
     *  - the array of local indexes might contain the same index multiple times as 
     *    a local entry might have halo copies on multiple processors.
     */
    HaloExchangePlan( 
        hmemo::HArray<IndexType> halo2GlobalIndexes,
        hmemo::HArray<IndexType> localIndexes,
        CommunicationPlan haloCommPlan, 
        CommunicationPlan localCommPlan,
        std::map<IndexType, IndexType> global2HaloMap );

    /** Constructor without global2Halo map, will be built */

    HaloExchangePlan( 
        hmemo::HArray<IndexType> halo2GlobalIndexes,
        hmemo::HArray<IndexType> localIndexes,
        CommunicationPlan haloCommPlan, 
        CommunicationPlan localCommPlan );

    /** Split up */

    void splitUp( hmemo::HArray<IndexType>& halo2GlobalIndexes,
                  hmemo::HArray<IndexType>& localIndexes,
                  CommunicationPlan& haloCommPlan,
                  CommunicationPlan& localCommPlan );

    /** Copy constructor. */

    HaloExchangePlan( const HaloExchangePlan& halo ) = default;

    /** Use default move contructor. */

    HaloExchangePlan( HaloExchangePlan&& halo ) = default;

    virtual ~HaloExchangePlan();

    /** Clear the halo plan, but do not free the allocated data */

    void clear();

    /** Purge, same as clear but also frees the allocated data */

    void purge();

    /** Swap member variables with another halo plan. */

    void swap( HaloExchangePlan& other );

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

    /** Provide updateHalo with a new haloArray, just for convenience 
     *
     *  This function with result argument is not recommended for multiple calls
     *  where it is more efficient to reuse objects that have capacity like HArray.
     */
    template<typename ValueType>
    hmemo::HArray<ValueType> updateHaloF(
        const hmemo::HArray<ValueType>& sourceArray, 
        const Communicator& comm ) const;

    /**
     *  same as updateHalo but provides a temporary array that takes the send buffer
     *
     *  Optimized version that reuses an allocated array for the send buffer.
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
     * 
     *  \code
     *     haloPlan.update( haloArray, sourceArray, comm, tmpSend );  
     *     // split it up as follows:
     *     utilskernel::HArrayUtils::gather( tmpSend, sourceArray, haloPlan.getLocalIndexes() );
     *     ....  // set up some asynchronous stuff
     *     haloPlan.updateDirect[Async]( haloArray, tmpSend, comm );
     *  \endcode
     *
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

    template<typename ValueType>
    tasking::SyncToken* updateHaloDirectAsync( 
        hmemo::HArray<ValueType>& haloArray, 
        const hmemo::HArray<ValueType>& sendArray, 
        const Communicator& comm ) const;
    /** 
     *  This method takes the halo (non-local values) to update the source array.
     */
    template<typename ValueType>
    void updateByHalo( 
        hmemo::HArray<ValueType>& sourceArray, 
        const hmemo::HArray<ValueType>& haloArray, 
        common::BinaryOp op, const Communicator& comm ) const;

    /** Use default assignment operator. */

    HaloExchangePlan& operator=( const HaloExchangePlan& ) = default;

    /** Use default move operator. */

    HaloExchangePlan& operator=( HaloExchangePlan&& ) = default;

    /** Getter for communication plan used into/from hala */

    inline const CommunicationPlan& getHaloCommunicationPlan() const;

    /** Getter for local communication plan used for local part of distributed array */

    inline const CommunicationPlan& getLocalCommunicationPlan() const;

    /** Getter for local indexes referring to local elements of distributed array used in other halos */

    inline const hmemo::HArray<IndexType>& getLocalIndexes() const;

    /**
     * Getter for the array that contains the translation from halo indexes to global indexes.
     */
    inline const hmemo::HArray<IndexType>& getHalo2GlobalIndexes() const;

    /** Query the size for a halo to be allocated 
     *
     *  This is an abbrevation either for getHalo2GlobalIndexes().size() or getHaloCommPlan().totalQuantitiy()
     */
    inline IndexType getHaloSize() const;

    /** 
     *  Predicate that returns true if this processor does not have to participate in update communication.
     *
     *  @returns true iff halo and local communication plans are empty.
     *
     *  Be careful: getHaloSize() == 0 implies that this processor has no counterparts of
     *  values, but it might be possible that this partition has local values with halo copies on other processors.
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

    /** 
     *   @brief Query the halo index for a (required) global index
     *
     *   @param[in] globalIndex is the global index that should have a halo counter part
     *   @returns position i with mHalo2Global[i] == globalIndex, invalidIndex otherwise
     */
    inline IndexType global2Halo( const IndexType globalIndex ) const;

    /** This method translates global indexes into local indexes 
     *
     *  @param[in] globalIndexes is the array with global indexes that are translated 
     *  @param[out] haloIndexes is the array with corresponding halo indexes.
     *
     *  The array must only contain required indexes used to build the halo.
     */
    void global2HaloV( hmemo::HArray<IndexType>& haloIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const;

    /** Constructor function as static version for more convenient logging */

    static HaloExchangePlan haloExchangePlan( 
        const class Distribution& distribution,
        const hmemo::HArray<IndexType>& globalIndexes, 
        const bool elimDouble );

private:

    // Indexes for required values and values to provide are stored in HArrays
    // so they might be used in different contexts, especially also on GPU

    hmemo::HArray<IndexType> mHalo2GlobalIndexes;   // non-local global indexes of the halo entries
    hmemo::HArray<IndexType> mLocalIndexes;         // local indexes of entries that have halo counterparts

    CommunicationPlan mHaloCommPlan;    // plan to communicate contiguous sections of halo
    CommunicationPlan mLocalCommPlan;   // is the transpose of mHaloCommPlan

    // This is an additional data structure to map the required global to halo indexes
    // mGlobal2Halo[mHalo2GlobalIndexes[i]] = i, 0 <= i < haloSize()

    std::map<IndexType, IndexType> mGlobal2Halo;  //!< inverse to mHalo2Global

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Help routine to eliminate double and local required indexes. */

    static void buildMap( 
        std::map<IndexType, IndexType>& globalMap,
        hmemo::HArray<PartitionId>& owners,
        const hmemo::HArray<IndexType>& globalIndexes,
        const IndexType NP );

    static void updateMap(
        std::map<IndexType,IndexType>& globalMap,
        const hmemo::HArray<IndexType>& globalIndexes );
};

const CommunicationPlan& HaloExchangePlan::getHaloCommunicationPlan() const
{
    return mHaloCommPlan;
}

const CommunicationPlan& HaloExchangePlan::getLocalCommunicationPlan() const
{
    return mLocalCommPlan;
}

const hmemo::HArray<IndexType>& HaloExchangePlan::getLocalIndexes() const
{
    return mLocalIndexes;
}

const hmemo::HArray<IndexType>& HaloExchangePlan::getHalo2GlobalIndexes() const
{
    return mHalo2GlobalIndexes;
}

IndexType HaloExchangePlan::global2Halo( const IndexType globalIndex ) const
{
    const std::map<IndexType, IndexType>::const_iterator elem = mGlobal2Halo.find( globalIndex );

    if ( elem == mGlobal2Halo.end() )
    {
        return invalidIndex;
    }

    return elem->second;
}

IndexType HaloExchangePlan::getHaloSize() const
{
    // Note: is also same as mHaloCommPlan.totalQuantity();

    return mHalo2GlobalIndexes.size();
}

bool HaloExchangePlan::isEmpty() const
{
    return ( mHaloCommPlan.totalQuantity() == 0 ) && ( mLocalCommPlan.totalQuantity() == 0 );
}

/** Build a halo exchange plan by an array of (unsorted) global indexes 
 *
 *  @param[in]  globalIndexes are global indexes for required values from other processors
 *  @param[in]  distribution is the mapping to find owners and map global to local indexes 
 *  @param[in]  elimDouble if true double entries in requiredIndexes are recognized
 *  @returns    the HaloExchangePlan object that contains (sorted) required and provides indexes and exchange plans
 *
 *  Notes:
 *
 *   - double entries in required indexes increase communication volume if elimDouble is false.
 *   - requiredIndexes usually does not contain global indexes that are local to this processor.
 *
 *  The halo2GlobalIndexes array of the plan is different from the requiredIndexes array. The first one
 *  is sorted by the owning processors while the entries in the latter one might be of any order.
 *  
 *  \code
 *      auto dist = blockDistribution( N );
 *      auto requiredIndexes = HArray<IndexType>( { ... } );  // global indexes of required data
 *      auto plan = haloExchangePlan( *dist, requiredIndexes );
 *      ...
 *      HArray<ValueType> localArray( dist->getLocalSize(), ... );   // distributed globalArray
 *      HArray<ValueType> haloArray( plan.getHaloSize() );           // halo part
 *      plan.updateHalo( haloArray, localArray );
 *      // for a halo index k haloArray[ k ] contains globalArray[ plan.halo2Global(k) ]
 *      // for a global index g haloArray[ plan.global2Halo(g) ] contains globalArray[g]
 *  \endcode
 */
inline HaloExchangePlan haloExchangePlan( 
    const class Distribution& distribution,
    const hmemo::HArray<IndexType>& globalIndexes, 
    const bool elimDouble = true )
{
    return HaloExchangePlan::haloExchangePlan( distribution, globalIndexes, elimDouble );
}

} /* end namespace dmemo */

} /* end namespace scai */
