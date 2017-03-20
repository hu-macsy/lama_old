/**
 * @file GeneralDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief GeneralDistribution.hpp
 * @author brandes
 * @date 25.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/utilskernel/LArray.hpp>
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace dmemo
{

/** A general distribution allows to map a global range of values
 to the partitions of a communicator completely arbitrarily.

 Each partition has the information which values it holds.
 A partition has no information where to find values
 not owned by itself.
 */

class COMMON_DLL_IMPORTEXPORT GeneralDistribution: public Distribution
{
public:

    /** Constructor of a general distribution where each processor knows it indexes.
     *
     *  \param globalSize is the size of the distributed range
     *  \param myGlobalIndexes contains all indexes of range owned by this partition
     *  \param communicator partitions on which the range is distributed.
     *
     *  Important: each global index from 0 to globalSize-1 must appear exactly once in
     *  the vector myGlobalIndexes on one partition.
     */
    GeneralDistribution(
        const IndexType globalSize,
        const hmemo::HArray<IndexType>& myGlobalIndexes,
        const CommunicatorPtr communicator );

    /** This constructor creates a general distribution by an array containing the owner for each element
     *
     *  @param[in] owners array with ower for each element, 0 <= owners[i] < communicator->size()
     *  @param[in] communicator that specifies the processor array for distribution
     *
     *  // Note: owners must only be valid on host processor
     */

    GeneralDistribution(
        const hmemo::HArray<PartitionId>& owners,
        const CommunicatorPtr communicator );

    explicit GeneralDistribution( const Distribution& other );

    /** Reimplment the default copy constructor */

    GeneralDistribution( const GeneralDistribution& other );

    virtual ~GeneralDistribution();

    /** Implementation of pure method Distribution::isLocal */

    virtual bool isLocal( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getLocalSize */

    virtual IndexType getLocalSize() const;

    /** Implementation of pure method Distribution::local2global */

    virtual IndexType local2global( const IndexType localIndex ) const;

    /** Implementation of pure method Distribution::global2local */

    virtual IndexType global2local( const IndexType globalIndex ) const;

    /** Implementation of pure function Distribution::getBlockDistributionSize.
     *
     *  Each processor must have a contiguous part of indexes and their order
     *  must match the rank order of processors.
     */
    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override the default implementation of Distribution::allOwners */

    virtual void allOwners( hmemo::HArray<PartitionId>& owners, const PartitionId root ) const;

    /** Override Distribution::getOwnedIndexes */

    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** Implementation of pure method Distribution::enableAnyAddressing */

    virtual void enableAnyAddressing() const;

    /** Implementation of pure method Distribution::getAnyLocalSize */

    virtual IndexType getAnyLocalSize( const PartitionId partition ) const;

    /** Implementation of pure method Distribution::getAnyOwner */

    virtual PartitionId getAnyOwner( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getAnyLocalIndex */

    virtual IndexType getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const;

    /** Implementation of pure method Distribution::getAnyGlobalIndex */

    virtual IndexType getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const;

    /** This method returns the array that contains for this processors all owned indexes. */

    inline const hmemo::HArray<IndexType>& getMyIndexes() const;

    /** Implementation of pure method Distribution::getKind */

    virtual inline const char* getKind() const;

protected:

    static const char theCreateValue[];

    /** This constructor might be called for derived classes that fill mGlobal2Local and mLocal2Global themselves. */

    GeneralDistribution( const IndexType globalSize, const CommunicatorPtr communicator );

    utilskernel::LArray<IndexType> mLocal2Global;   //!< for each local index its global index, entries are sorted
 
    // the following arrays will only be available if enableAnyAddressing has been called
    // Note: if set the array mGlobal2Local is no more needed

    mutable utilskernel::LArray<PartitionId> mAllOwners;
    mutable utilskernel::LArray<IndexType> mAllLocalOffsets;     // local size on each partition
    mutable utilskernel::LArray<IndexType> mAllLocal2Global;     // sorts elements into buckets 
    mutable utilskernel::LArray<IndexType> mAllGlobal2Local;     // sorts elements into buckets 

    // Example
    // index       0    1    2    3   4    5    6    7   8   9   10   11   12 
    // mOwners:    0    1    2    0   2    0    1    0   0   1    1    2    2 
    // Offsets:    0                       5                 9                    13
    // perm   :    0    3    5    7   8    1    6    9  10   2    4   11   12     local2global
    // perm'  :    0    5    9    1  10    2    6    3   4   7    8   11   12     global2local
    // 
    // Note: perm is identity iff we have a block distribution

private:

    GeneralDistribution();

    GeneralDistribution& operator=( const GeneralDistribution& other );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

typedef common::shared_ptr<GeneralDistribution> GeneralDistributionPtr;

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

const hmemo::HArray<IndexType>& GeneralDistribution::getMyIndexes() const
{
    return mLocal2Global;
}

const char* GeneralDistribution::getKind() const
{
    return theCreateValue;
}

} /* end namespace dmemo */

} /* end namespace scai */
