/**
 * @file include/scai/testsupport/uniquePathComm.hpp
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
 * @brief Unique path generation independent/shared among nodes in a communicator.
 * @author Andreas Longva
 * @date 21.11.2017
 */
#pragma once

#include <scai/dmemo/Communicator.hpp>
#include <scai/testsupport/uniquePath.hpp>

#include <string>
#include <sstream>

namespace scai
{

namespace testsupport
{

    /**
     * Generate a unique path which is unique also among the nodes in the communicator.
     *
     * See uniquePath() for the role of `dir` and `namePrefix`.
     */
    inline std::string uniquePathPerNode(const std::string & dir,
                                  const scai::dmemo::Communicator & comm,
                                  const std::string & namePrefix = "")
    {
        std::stringstream path;
        path << uniquePath(dir, namePrefix);
        path << '_' << comm.getSize() << '_' << comm.getRank();
        return path.str();
    }

    /**
     * Generate a unique path which is shared among the nodes in the communicator.
     *
     * See uniquePath() for the role of `dir` and `namePrefix`.
     */
    inline std::string uniquePathSharedAmongNodes(const std::string & dir,
                                           const scai::dmemo::Communicator & comm,
                                           const std::string & namePrefix = "")
    {
        // TODO: Make exception-safe with respect to MPI
        std::string path;
        if (comm.getRank() == 0)
        {
            path = ::scai::testsupport::uniquePath(dir, namePrefix);
        }
        comm.bcast(path, 0);
        return path;
    }

} // namespace testsupport

} // namespace scai
