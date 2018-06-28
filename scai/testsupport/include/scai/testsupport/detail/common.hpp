/**
 * @file include/scai/testsupport/detail/common.hpp
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
 * @brief Common detail functionality for testsupport code.
 * @author Andreas Longva
 * @date 09.11.2017
 */
#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

// Need the following macro to retrieve the BOOST_TEST_MODULE test name
#define LAMATEST_STRINGIFY2(x) #x
#define LAMATEST_STRINGIFY(x) LAMATEST_STRINGIFY2(x)

namespace scai
{

namespace testsupport
{

namespace detail
{

inline std::string boostTestModuleName()
{
    return std::string(LAMATEST_STRINGIFY(BOOST_TEST_MODULE));
}

inline std::pair<std::string, std::string> getKeyValueFromArg(const std::string & arg)
{
    const auto equalsPos = std::find(arg.begin(), arg.end(), '=');
    std::string key, value;
    std::copy(arg.begin(), equalsPos, std::back_inserter(key));
    std::copy(equalsPos + 1, arg.end(), std::back_inserter(value));
    return std::make_pair(std::move(key), std::move(value));
}

inline bool argIsSink(const std::string & arg, const std::string & sinkName)
{
    std::string key;
    std::tie(key, std::ignore) = getKeyValueFromArg(arg);
    return key == sinkName;
}

inline std::vector<char> intoNullTerminatedChars(const std::string & str)
{
    std::vector<char> chars;
    std::copy(str.begin(), str.end(), std::back_inserter(chars));
    chars.push_back('\0');
    return chars;
}

struct ArgParseResult
{
    std::string tempDir;
    std::vector<std::vector<char>> args;
};

inline ArgParseResult parseAndRebuildArgs(int argc, char ** argv, const std::string & testSuiteName)
{
    static const std::string REPORT_SINK_ARG = "--report_sink";
    static const std::string LOG_SINK_ARG = "--log_sink";
    static const std::string OUTPUT_DIR_ARG = "--output_dir";
    static const std::string TEMP_DIR_ARG = "--temp_dir";

    ArgParseResult result;

    const std::array<std::string, 4> forbiddenSinks {
        { REPORT_SINK_ARG, LOG_SINK_ARG, std::string("-e"), std::string("-k") }
    };

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

        // Do not include non-Boost options in new list of arguments,
        // since Boost does not recognize them and will complain
        if (key == OUTPUT_DIR_ARG)
        {
            outputDir = value;
        }
        else if (key == TEMP_DIR_ARG)
        {
            result.tempDir = value;
        }
        else
        {
            result.args.push_back(intoNullTerminatedChars(arg));
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
        result.args.push_back(intoNullTerminatedChars(report_sink));
        result.args.push_back(intoNullTerminatedChars(log_sink));
    }

    return result;
}

} // namespace detail

} // namespace testsupport

} //namespace scai
