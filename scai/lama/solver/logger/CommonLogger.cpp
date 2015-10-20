/**
 * @file CommonLogger.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the CommonLogger class
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/logger/CommonLogger.hpp>

namespace scai
{

namespace lama
{

CommonLogger::CommonLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    bool ignoreRank )
    : SolverLogger( id, level, writeBehaviour, ignoreRank )
{
}

CommonLogger::CommonLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    common::shared_ptr<Timer> timer,
    bool ignoreRank )
    : SolverLogger( id, level, writeBehaviour, timer, ignoreRank )
{
}

CommonLogger::CommonLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    const std::string& logFileName,
    common::shared_ptr<Timer> timer,
    bool ignoreRank )
    : SolverLogger( id, level, writeBehaviour, logFileName, timer, ignoreRank )
{
}

CommonLogger::~CommonLogger()
{
}

std::string CommonLogger::createPrefix()
{
    return mId;
}

} /* end namespace lama */

} /* end namespace scai */
