/**
 * @file include/scai/testsupport/uniquePath.hpp
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
 * @brief Facilities for unique path generation for tests (temporary path generation).
 * @author Andreas Longva
 * @date 21.11.2017
 */
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
