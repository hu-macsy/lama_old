/**
 * @file GeneralDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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

// std
#include <vector>
#include <map>

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

    /** This constructor creates a general distribution by a mapping of rows to partition ids.
     *
     *  @param[in] row2Partition has globalSize entries, row2Partition[i] specifies owner of i
     *  @param[in] globalSize is the number of distributed elements
     *  @param[in] communicator that specifies the processor array for distribution
     *
     *  Note:  0 <= row2Partion[i] < communicator->size() for 0 <= i < globalSize
     *  Note:  Only the master process has to provide the mapping row2Partition.
     */

    /*
    GeneralDistribution(
        const std::vector<IndexType>& row2Partition,
        const IndexType globalSize,
        const CommunicatorPtr communicator );
    */

    /** This constructor creates a general distribution by an array containing the owner for each element
     *
     *  @param[in] owners, with 0 <= owners[i] < communicator->size()
     *  @param[in] communicator that specifies the processor array for distribution
     *
     *  // Note: owners must only be valid on host processor
     */

    GeneralDistribution(
        const hmemo::HArray<IndexType>& owners,
        const CommunicatorPtr communicator );

    explicit GeneralDistribution( const Distribution& other );

//    GeneralDistribution(const GeneralDistribution& other);

    virtual ~GeneralDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual IndexType getLocalSize() const;

    // virtual std::vector<IndexType>& getLocalRows();

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override the default implementation of Distribution::allOwners */

    virtual void allOwners( hmemo::HArray<PartitionId>& owners, const PartitionId root ) const;

    /** Override Distribution::getOwnedIndexes */

    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** This method returns the array that contains for this processors all owned indexes. */

    inline const hmemo::HArray<IndexType>& getMyIndexes();

    virtual const char* getKind() const
    {
        return theCreateValue;
    }

protected:

    static const char theCreateValue[];

    GeneralDistribution( const IndexType globalSize, const CommunicatorPtr communicator );

    typedef std::map<IndexType, IndexType> Global2LocalMapType;

    Global2LocalMapType mGlobal2Local;

    utilskernel::LArray<IndexType> mLocal2Global;

private:

    GeneralDistribution();

    GeneralDistribution& operator=( const GeneralDistribution& other );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

typedef common::shared_ptr<GeneralDistribution> GeneralDistributionPtr;

const hmemo::HArray<IndexType>& GeneralDistribution::getMyIndexes()
{
    return mLocal2Global;
}

} /* end namespace dmemo */

} /* end namespace scai */
