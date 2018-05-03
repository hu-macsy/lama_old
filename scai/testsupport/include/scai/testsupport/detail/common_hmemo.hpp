/**
 * @file include/scai/testsupport/detail/common_hmemo.hpp
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
 * @brief Common detail functionality for testsupport code that rely on hmemo.
 * @author Andreas Longva
 * @date 09.11.2017
 */
#pragma once

#include <string>
#include <scai/hmemo/Context.hpp>

namespace scai
{

namespace testsupport
{

namespace detail
{

/**
 * Find a suitable name for a test suite which describes
 * the context used. Note that the naming is chosen such that
 * related environments are collected when sorted.
 */
std::string adaptTestSuiteNameToEnv(const std::string & name, const hmemo::Context & context)
{
    std::string prefix;
    switch (context.getType())
    {
        case common::ContextType::Host:
            prefix = "~Host ";
            break;
        case common::ContextType::CUDA:
            prefix = "~CUDA ";
            break;
        default:
            std::cerr << "Unsupported context type. Can not create appropriate test suite name." << std::endl;
            throw std::runtime_error("Unsupported context type.");
    }

    return prefix + name;
}

std::string suiteNameForFile(const std::string & name, const hmemo::Context & context)
{
    // TODO: Context
    std::stringstream filename;
    filename << name;

    switch (context.getType())
    {
        case common::ContextType::Host:
            filename << "_host";
            break;
        case common::ContextType::CUDA:
            filename << "_cuda";
            break;
        default:
            std::cerr << "Unsupported context type. Can not create appropriate test filename." << std::endl;
            throw std::runtime_error("Unsupported context type.");
    }

    return filename.str();
}

} // namespace detail

} // namespace testsupport

} // namespace scai
