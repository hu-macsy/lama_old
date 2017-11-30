#pragma once

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

#include <scai/common/Settings.hpp>

#include <scai/testsupport/detail/common.hpp>

namespace scai
{

namespace testsupport
{

int commonTestMain( int argc, char* argv[] )
{
    using scai::testsupport::detail::rebuildArgs;

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    const auto boostTestModuleName = std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
    const auto testSuiteName = boostTestModuleName;

    // Building args as a vector<vector<char>> ensures that lifetime of modified args is bounded by main() call
    auto newArgs = rebuildArgs(argc, argv, testSuiteName);
    std::vector<char *> charPointers;
    for (auto & arg : newArgs)
    {
        charPointers.push_back(arg.data());
    }

    return boost::unit_test::unit_test_main(&init_unit_test, charPointers.size(), charPointers.data());
}

} // namespace testsupport

} // namespace scai
