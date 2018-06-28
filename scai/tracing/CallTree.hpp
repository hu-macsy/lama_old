/**
 * @file CallTree.hpp
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
 * @brief Definition of static class that provides routines to generate a Call Tree.
 * @author Thomas Brandes
 * @date 03.07.2011
 */

#pragma once

// local library
#include <scai/tracing/RegionEntry.hpp>
#include <scai/tracing/Counters.hpp>

// internal scai libraries
#include <scai/logging.hpp>

// std
#include <vector>

namespace scai
{

namespace tracing
{

class CallTree
{
public:

    static void enter( const int region_id, RegionEntry& region, const CounterArray& startValues );

    static void leave( const int region_id, const RegionEntry& region, const CounterArray& stopValues );

    static void finish();

private:

    static std::vector<int> theCallStack;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace tracing */

} /* end namespace scai */
