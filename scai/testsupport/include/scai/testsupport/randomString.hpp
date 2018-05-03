/**
 * @file include/scai/testsupport/randomString.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Functionality for random string generation.
 * @author Andreas Longva
 * @date 21.11.2017
 */
#pragma once

#include <random>
#include <string>
#include <cassert>

namespace scai
{

namespace testsupport
{

template <typename RandomEngine>
std::string randomAlphaNumericString(RandomEngine & engine, size_t len)
{
    const static std::string allowed_chars =
        "01234567890"
        "abcdefghjiklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    // Need to use unsigned int instead of size_t because it is not guaranteed
    // that a distribution for size_t exists
    std::uniform_int_distribution<unsigned int> distribution(
        static_cast<unsigned int>(0),
        static_cast<unsigned int>(allowed_chars.size() - 1)
    );

    std::string s;
    s.reserve(len);
    for (size_t i = 0; i < len; ++i)
    {
        const auto random_index = distribution(engine);
        assert(random_index < allowed_chars.size());
        const auto random_char = allowed_chars[static_cast<size_t>(random_index)];
        s.push_back(random_char);
    }

    return s;
}

} // namespace testsupport

} // namespace scai
