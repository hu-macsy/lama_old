/**
 * @file include/scai/testsupport/dmemoTestMain.hpp
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
 * @brief Test binary main() function for dmemo-dependent test binaries.
 * @author Andreas Longva
 * @date 09.11.2017
 */
#pragma once

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/dmemo.hpp>
#include <scai/testsupport/detail/common.hpp>
#include <scai/testsupport/detail/common_hmemo.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

#include <iostream>
#include <sstream>


namespace scai
{

namespace testsupport
{

namespace detail
{

/**
 * Find a suitable name for a test suite which describes
 * the context/communicator used. Note that the naming is chosen such that
 * related environments are collected when sorted, rather than
 * related tests. For example, all MPI tests which have a certain number of processors
 * should appear clustered when test suite names are sorted.
 */
std::string adaptTestSuiteNameToEnv(const std::string & name,
                                    const scai::hmemo::Context & context,
                                    const scai::dmemo::Communicator & comm)
{
    // TODO: Context
    std::stringstream newTestName;
    if (comm.getType() == scai::dmemo::_Communicator::MPI)
    {
        newTestName << "~~MPI " << comm.getSize() << ":" << comm.getRank() << " ";
    }

    newTestName << adaptTestSuiteNameToEnv(name, context);
    return newTestName.str();
}

std::string suiteNameForFile(const std::string & suiteName,
                             const scai::hmemo::Context & context,
                             const scai::dmemo::Communicator & comm)
{
    std::stringstream filename;
    filename << suiteNameForFile(suiteName, context);
    if (comm.getType() == scai::dmemo::_Communicator::MPI)
    {
        filename << "_mpi_" << comm.getSize() << "_" << comm.getRank();
    }
    return filename.str();
}

bool dmemo_test_init()
{
    int nThreads;

    if ( scai::common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    scai::hmemo::ContextPtr ctx;
    try
    {
        ctx = scai::hmemo::Context::getContextPtr();
    }
    catch ( scai::common::Exception& ex )
    {
        std::cerr << "Could not get context for test: " << ex.what() << std::endl;
        return false;
    }

    scai::dmemo::CommunicatorPtr comm;
    try
    {
        comm = scai::dmemo::Communicator::getCommunicatorPtr();
    }
    catch ( scai::common::Exception& ex )
    {
        std::cerr << "Could not get the default communicator: " << ex.what() << std::endl;
        return false;
    }

    auto & master_suite = boost::unit_test::framework::master_test_suite();
    const auto suiteName = boostTestModuleName();
    const auto newTestName = adaptTestSuiteNameToEnv(suiteName, *ctx, *comm);
    master_suite.p_name.value = newTestName;

    return true;
}



} // namespace detail

int dmemoTestMain( int argc, char* argv[] )
{
    using scai::testsupport::detail::suiteNameForFile;
    using scai::testsupport::detail::parseAndRebuildArgs;
    using scai::testsupport::detail::dmemo_test_init;

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    const auto ctx = scai::hmemo::Context::getContextPtr();
    const auto comm = scai::dmemo::Communicator::getCommunicatorPtr();
    const auto boostTestModuleName = std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
    const auto testSuiteName = suiteNameForFile(boostTestModuleName, *ctx, *comm);

    // Building args as a vector<vector<char>> ensures that lifetime of modified args is bounded by main() call
    auto parseResult = parseAndRebuildArgs(argc, argv, testSuiteName);
    std::vector<char *> charPointers;
    for (auto & arg : parseResult.args)
    {
        charPointers.push_back(arg.data());
    }

    GlobalTempDir::setPathOrDefault(parseResult.tempDir);

    return boost::unit_test::unit_test_main( &dmemo_test_init, charPointers.size(), charPointers.data() );
}

} // namespace testsupport

} // namespace scai
