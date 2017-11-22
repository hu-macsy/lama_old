#pragma once

#include <string>
#include <random>
#include <chrono>
#include <type_traits>

namespace scai
{

namespace testsupport
{

/**
 * Generate a path extremely likely to be unique.
 *
 * Generates a path in the given directory path which is extremely likely to be unique.
 * The pattern of the generated path is given by:
 *
 * dir/namePrefixRANDOM
 *
 * where RANDOM is a random sequence of alpha-numeric characters of unspecified length.
 *
 * Note: The path manipulation is extremely rudimentary and entirely string-based.
 */
std::string uniquePath(const std::string & dir, const std::string & namePrefix = "");

} // namespace testsupport

} // namespace scai
