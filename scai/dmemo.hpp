/**
 * @file dmemo.hpp
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
 * @brief General header file that includes the main header files of the dmemo subproject.
 * @author Thomas Brandes
 * @date 09.02.2016
 */

// We include only base classes, not derived Communicator or Distribution classes

#include <scai/dmemo/Communicator.hpp>

#include <scai/dmemo/HaloExchangePlan.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
