/**
 * @file VTInterface.hpp
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
 * @brief Interface between this Trace library and VampirTrace
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

namespace scai
{

namespace tracing
{

class RegionEntry;

/** Static class provding interface to Vampir */

class VTInterface
{
public:

    /** Type definition needed as region entry gets additional entry to save region id given by VampirTrace. */

    typedef unsigned int VTRegionId;

    /** This routine defines a new region within VampirTrace. */

    static void define( RegionEntry& );

    static void enter( const RegionEntry& region );

    static void leave( const RegionEntry& region );

    static void enable( bool flag );
};

} /* end namespace tracing */

} /* end namespace scai */
