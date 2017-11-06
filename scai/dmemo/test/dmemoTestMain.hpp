#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/dmemo.hpp>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <array>

// Need the following macro to retrieve the BOOST_TEST_MODULE test name
#define LAMATEST_STRINGIFY2(x) #x
#define LAMATEST_STRINGIFY(x) LAMATEST_STRINGIFY2(x)

namespace dmemoMainDetail
{

/**
 * Find a suitable name for a test suite which describes
 * the context/communicator used. Note that the naming is chosen such that
 * related environments are collected when sorted, rather than
 * related tests. For example, all MPI tests which have a certain number of processors
 * should appear clustered when test suite names are sorted.
 */
std::string adaptTestSuiteNameToEnv(const std::string & name, const scai::dmemo::Communicator & comm)
{
    // TODO: Context
    std::stringstream newTestName;
    if (comm.getType() == scai::dmemo::_Communicator::MPI)
    {
        newTestName << "_MPI (" << comm.getSize() << " procs, rank: " << comm.getRank() << ") ";
    }

    newTestName << name;
    return newTestName.str();
}

std::string suiteNameForFile(const std::string & suiteName,
                             const scai::dmemo::Communicator & comm)
{
    // TODO: Context
    std::stringstream filename;
    filename << suiteName;
    if (comm.getType() == scai::dmemo::_Communicator::MPI)
    {
        filename << "_mpi_" << comm.getSize() << "_" << comm.getRank();
    }
    return filename.str();
}

bool init_function()
{
    int nThreads;

    if ( scai::common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
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
    const auto suiteName = std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
    const auto newTestName = adaptTestSuiteNameToEnv(suiteName, *comm);
    master_suite.p_name.value = newTestName;

    return true;
}

std::pair<std::string, std::string> getKeyValueFromArg(const std::string & arg)
{
    const auto equalsPos = std::find(arg.begin(), arg.end(), '=');
    std::string key, value;
    std::copy(arg.begin(), equalsPos, std::back_inserter(key));
    std::copy(equalsPos + 1, arg.end(), std::back_inserter(value));
    return std::make_pair(std::move(key), std::move(value));
}

bool argIsSink(const std::string & arg, const std::string & sinkName)
{
    std::string key;
    std::tie(key, std::ignore) = getKeyValueFromArg(arg);
    return key == sinkName;
}

std::vector<char> intoNullTerminatedChars(const std::string & str)
{
    std::vector<char> chars;
    std::copy(str.begin(), str.end(), std::back_inserter(chars));
    chars.push_back('\0');
    return chars;
}

std::vector<std::vector<char>> rebuildArgs(int argc, char ** argv, const std::string & testSuiteName)
{
    static const std::string REPORT_SINK_ARG = "--report_sink";
    static const std::string LOG_SINK_ARG = "--log_sink";
    static const std::string OUTPUT_DIR_ARG = "--output_dir";

    const std::array<std::string, 4> forbiddenSinks {
        { REPORT_SINK_ARG, LOG_SINK_ARG, std::string("-e"), std::string("-k") }
    };

    std::vector<std::vector<char>> args;
    std::string outputDir;
    bool sinkIsPresent = false;
    for (int i = 0; i < argc; ++i)
    {
        const auto arg = std::string(argv[i]);
        std::string key, value;
        std::tie(key, value) = getKeyValueFromArg(arg);

        const auto argHasSink = std::any_of(forbiddenSinks.begin(), forbiddenSinks.end(),
                    [&key] (const std::string & sink) { return sink == key; });
        sinkIsPresent = sinkIsPresent || argHasSink;

        if (key == OUTPUT_DIR_ARG)
        {
            // Do not include in new list of arguments, since Boost does not recognize it and will complain
            outputDir = value;
        } else {
            args.push_back(intoNullTerminatedChars(arg));
        }
    }

    if (sinkIsPresent && !outputDir.empty())
    {
        std::cerr << "Report/log sinks are not allowed when "
                  << OUTPUT_DIR_ARG << " has been supplied." << std::endl;
        throw std::runtime_error("Unsupported command-line argument.");
    }

    if (!outputDir.empty())
    {
        const auto report_sink = REPORT_SINK_ARG + "=" + outputDir + "/" + testSuiteName + "_report.xml";
        const auto log_sink = LOG_SINK_ARG + "=" + outputDir + "/" + testSuiteName + "_log.xml";
        args.push_back(intoNullTerminatedChars(report_sink));
        args.push_back(intoNullTerminatedChars(log_sink));
    }

    return args;
}

}

int dmemoTestMain( int argc, char* argv[] )
{
    using namespace dmemoMainDetail;

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    const auto comm = scai::dmemo::Communicator::getCommunicatorPtr();
    const auto boostTestModuleName = std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
    const auto testSuiteName = suiteNameForFile(boostTestModuleName, *comm);

    // Building args as a vector<vector<char>> ensures that lifetime of modified args is bounded by main() call
    auto newArgs = rebuildArgs(argc, argv, testSuiteName);
    std::vector<char *> charPointers;
    for (auto & arg : newArgs)
    {
        charPointers.push_back(arg.data());
    }

    return boost::unit_test::unit_test_main( &init_function, charPointers.size(), charPointers.data() );
}
