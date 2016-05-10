/**
 * @file NoDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

    NoDistribution( const IndexType globalSize );

    virtual ~NoDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual IndexType getLocalSize() const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    void printDistributionVector( std::string name ) const;

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

    virtual const char* getKind() const
    {
        return theCreateValue;
    }

private:

    NoDistribution(); // no default constructor as global size is not available

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static const char theCreateValue[];
};

} /* end namespace dmemo */

} /* end namespace scai */
