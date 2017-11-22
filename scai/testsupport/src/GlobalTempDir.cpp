#include <scai/testsupport/GlobalTempDir.hpp>

#include <memory>
#include <string>

namespace scai
{

namespace testsupport
{

std::unique_ptr<std::string> GlobalTempDir::m_tempDirPath = std::unique_ptr<std::string>();

}

}
