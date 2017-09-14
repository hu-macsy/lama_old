/**
 * @file BenchmarkPrinter.cpp
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
 * @brief BenchmarkPrinter.cpp
 * @author Jiri Kraus
 * @date 13.05.2011
 */
#include <scai/benchmark/BenchmarkPrinter.hpp>

namespace bf
{

bool BenchmarkPrinter::doOutput = true;

void BenchmarkPrinter::setDoOutput( const bool doOutput )
{
    BenchmarkPrinter::doOutput = doOutput;
}

const char* const BenchmarkPrinter::message_begin = "%_BENCHMARK_FRAMEWORK_!MESSAGE_START_!";
const char* const BenchmarkPrinter::message_end = "%_BENCHMARK_FRAMEWORK_!MESSAGE_END_!";
const char* const BenchmarkPrinter::error_begin = "%_BENCHMARK_FRAMEWORK_!ERROR_START_!";
const char* const BenchmarkPrinter::error_end = "%_BENCHMARK_FRAMEWORK_!ERROR_END_!";
const char* const BenchmarkPrinter::warning_begin = "%_BENCHMARK_FRAMEWORK_!WARNING_START_!";
const char* const BenchmarkPrinter::warning_end = "%_BENCHMARK_FRAMEWORK_!WARNING_END_!";

}
