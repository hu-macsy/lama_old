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
    std::string uniquePathPerNode(const std::string & dir,
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
    std::string uniquePathSharedAmongNodes(const std::string & dir,
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
