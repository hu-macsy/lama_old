/**
 * @file include/scai/testsupport/commonTestMain.hpp
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
 * @brief Test binary main() function for non-heterogeneous, non-distributed tests.
 * @author Andreas Longva
 * @date 09.11.2017
 */
#pragma once

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

#include <scai/common/Settings.hpp>

#include <scai/testsupport/detail/common.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

namespace scai
{

namespace testsupport
{

int commonTestMain( int argc, char* argv[] )
{
    using scai::testsupport::detail::parseAndRebuildArgs;

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    const auto boostTestModuleName = std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
    const auto testSuiteName = boostTestModuleName;

    // Building args as a vector<vector<char>> ensures that lifetime of modified args is bounded by main() call
    auto parseResult = parseAndRebuildArgs(argc, argv, testSuiteName);
    std::vector<char *> charPointers;
    for (auto & arg : parseResult.args)
    {
        charPointers.push_back(arg.data());
    }

    GlobalTempDir::setPathOrDefault(parseResult.tempDir);

    return boost::unit_test::unit_test_main(&init_unit_test, charPointers.size(), charPointers.data());
}

} // namespace testsupport

} // namespace scai
