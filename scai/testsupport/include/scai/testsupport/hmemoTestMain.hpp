#pragma once

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/hmemo.hpp>
#include <scai/testsupport/detail/common.hpp>
#include <scai/testsupport/detail/common_hmemo.hpp>

#include <iostream>
#include <sstream>
#include <exception>

namespace scai
{

namespace testsupport
{

namespace detail
{

bool hmemo_test_init()
{
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

    int nThreads;
    if ( scai::common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    auto & master_suite = boost::unit_test::framework::master_test_suite();
    const auto suiteName = boostTestModuleName();
    const auto newTestName = adaptTestSuiteNameToEnv(suiteName, *ctx);
    master_suite.p_name.value = newTestName;

    return true;
}



} // namespace detail

int hmemoTestMain( int argc, char* argv[] )
{
    using scai::testsupport::detail::boostTestModuleName;
    using scai::testsupport::detail::suiteNameForFile;
    using scai::testsupport::detail::rebuildArgs;
    using scai::testsupport::detail::hmemo_test_init;

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    const auto ctx = scai::hmemo::Context::getContextPtr();
    const auto testSuiteName = suiteNameForFile(boostTestModuleName(), *ctx);

    // Building args as a vector<vector<char>> ensures that lifetime of modified args is bounded by main() call
    auto newArgs = rebuildArgs(argc, argv, testSuiteName);
    std::vector<char *> charPointers;
    for (auto & arg : newArgs)
    {
        charPointers.push_back(arg.data());
    }

    return boost::unit_test::unit_test_main( &hmemo_test_init, charPointers.size(), charPointers.data() );
}

} // namespace testsupport

} // namespace scai