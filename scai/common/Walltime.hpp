/**
 * @file common/Walltime.hpp
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
 * @brief Definition of static class that provides high precision walltime
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// std
#include <stdint.h>

namespace scai
{

namespace common
{

/** Make a common defintion for 8-byte integer as uint64_t can be implementation specific. */

typedef uint64_t INTEGER_8;

/**
 * @brief A simple static class that delivers walltime (used for logging and tracing)
 */
class COMMON_DLL_IMPORTEXPORT Walltime
{
public:

    /** Get the current walltime.
     *
     *  @return current walltime in ticks
     */
    static INTEGER_8 timestamp();

    /** Number of ticks per second, so timestamp() / timerate() gives time in seconds. */

    static INTEGER_8 timerate();

    /** Get the current walltime.
     *
     *  @return current walltime in seconds
     */
    static double get();

    /**
     *  sleep routine for milliseconds
     */

    static void sleep( unsigned int milliseconds );

private:

    /** Private constructor for a static class. */

    Walltime();
};

} /* end namespace common */

} /* end namespace scai */
