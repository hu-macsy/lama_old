/**
 * @file Walltime.hpp
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
 * @brief Class that gives back walltime
 * @author Thomas Brandes
 * @date 25.04.2013
 */

#pragma once

namespace lama
{

/**
 * @brief A simple static class that delivers walltime (used for logging and tracing)
 */
class Walltime
{
public:

    /** Get the current walltime.
     *
     *  @return current walltime in seconds
     */
    static double get();

private:

    /** Private constructor for a static class. */

    Walltime();
};

} // namespace lama
