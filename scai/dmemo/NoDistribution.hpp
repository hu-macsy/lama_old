/**
 * @file NoDistribution.hpp
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
 * @brief NoDistribution.hpp
 * @author brandes
 * @date 14.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/logging.hpp>

namespace scai
{

namespace dmemo
{

/** Distribution class that stands for a replicated distribution.
 *
 *  With this distribution an object has an incarnation on each
 *  processor.
 *
 *  Usually, methods should take care of consistency among
 *  all processors, i.e. writes and update operations must be
 *  done on all partitions. But a replicated object can also be used
 *  like a private incarnation on each processor.
 */
class COMMON_DLL_IMPORTEXPORT NoDistribution:

    public Distribution,
    private Distribution::Register<NoDistribution>

{
public:

    /** Constructor of NoDistribution 
     *
     *  @param[in] globalSize is the global size of the distributed object
     *  @param[in] comm       specifies the set of processors that have an incarnation 
     */
    NoDistribution( const IndexType globalSize, CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

    virtual ~NoDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual IndexType getLocalSize() const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    /** Implementation of pure function Distribution::getBlockDistributionSize, here same as getLocalSize */

    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static inline std::string createValue();

    virtual inline const char* getKind() const;

    static const char* getId();

    /** Implementation of pure method Distribution::hasAnyAddressing */
    virtual bool hasAnyAddressing() const;

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

private:

    NoDistribution(); // no default constructor as global size is not available

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

/* ------------------------------------------------------------------------ */
/*  Implementation of inline methods                                        */
/* ------------------------------------------------------------------------ */

const char* NoDistribution::getKind() const
{
    return getId();
}

std::string NoDistribution::createValue()
{
    return getId();
}

} /* end namespace dmemo */

} /* end namespace scai */
