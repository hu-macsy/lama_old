/**
 * @file BlockPartitioning.hpp
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
 * @brief BlockPartitioning.hpp
 * @author Thomas Brandes
 * @date 18.08.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/partitioning/Partitioning.hpp>

// std
#include <vector>

namespace scai
{

namespace partitioning
{

/** This block partitioning creates just a block distribution where the size of the local part depends on its weight.
 *
 */

class COMMON_DLL_IMPORTEXPORT BlockPartitioning:

    public Partitioning,
    private Partitioning::Register<BlockPartitioning>

{
public:

    /** Constructor */

    BlockPartitioning();

    /** Destructor */

    virtual ~BlockPartitioning();

    /** Implementation of pure method Partitioning::rectangularPartitioning */

    virtual void rectangularPartitioning( hmemo::HArray<PartitionId>& rowMapping,
                                          hmemo::HArray<PartitionId>& colMapping,
                                          const lama::_Matrix& matrix,
                                          const hmemo::HArray<float>& processorWeights ) const;

    /** Implementation of pure method Partitioning::squarePartitioning */

    virtual void squarePartitioning( hmemo::HArray<PartitionId>& newLocalOwners, 
                                     const lama::_Matrix& matrix, 
                                     const hmemo::HArray<float>& processorWeights ) const;

    /** Override of Partitioning::rectangularRedistribute */

    virtual void rectangularRedistribute( lama::_Matrix& matrix, const float weight ) const;

    /** Override Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /** Static method required for create to use in Partitioning::Register */

    static PartitioningPtr create();

    /** Static method required for Partitioning::Register */

    static std::string createValue();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    static const char theCreateValue[];
};

} /* end namespace partitioning */

} /* end namespace scai */
