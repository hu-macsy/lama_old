/**
 * @file CommonLogger.cpp
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
 * @brief Contains the implementation of the CommonLogger class
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/logger/CommonLogger.hpp>

namespace scai
{

namespace solver
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
    std::shared_ptr<Timer> timer,
    bool ignoreRank )
    : SolverLogger( id, level, writeBehaviour, timer, ignoreRank )
{
}

CommonLogger::CommonLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    const std::string& logFileName,
    std::shared_ptr<Timer> timer,
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

} /* end namespace solver */

} /* end namespace scai */
