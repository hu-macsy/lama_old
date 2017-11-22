#pragma once

#include <scai/dmemo/Communicator.hpp>
#include <scai/testsupport/uniquePath.hpp>

#include <string>
#include <sstream>

namespace scai
{

namespace testsupport
{

    std::string uniquePathPerNode(const std::string & unique_dir,
                                  const scai::dmemo::Communicator & comm,
                                  const std::string & namePrefix = "")
    {
        std::stringstream path;
        path << uniquePath(unique_dir, namePrefix);
        path << '_' << comm.getSize() << '_' << comm.getRank();
        return path.str();
    }

    std::string uniquePathSharedAmongNodes(const std::string & unique_dir,
                                           const scai::dmemo::Communicator & comm,
                                           const std::string & namePrefix = "")
    {
        // TODO: Make exception-safe
        std::string path;
        if (comm.getRank() == 0)
        {
            path = ::scai::testsupport::uniquePath(unique_dir, namePrefix);
        }
        comm.bcast(path, 0);
        return path;
    }

} // namespace testsupport

} // namespace scai
