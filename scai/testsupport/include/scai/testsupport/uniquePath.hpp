#pragma once

#include <string>
#include <random>
#include <chrono>
#include <type_traits>

namespace scai
{

namespace testsupport
{

std::string uniquePath(const std::string & dir, const std::string & namePrefix = "");

} // namespace testsupport

} // namespace scai
