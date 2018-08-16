/**
 * @file AbstractLoggerCreator.hpp
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
 * @brief Abstract class for logger creator.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// local library
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace logging
{

/** Abstract class from which classes must derive that can create loggers. */

class AbstractLoggerCreator
{

public:
    /**
     * @brief Destructor needed due to virtual functions
     */
    virtual ~AbstractLoggerCreator()
    {
    }

    /**
     * @brief get the root logger
     */

    virtual class Logger& getRoot() const = 0;

    /**
     * @brief Function to create a new instance at a given level.
     * @param[in] name is identification of logger, must not contain any dots
     * @param[in,out] parent is the parent logger, new logger is set as son
     */

    virtual class Logger* create( const std::string& name, class Logger* parent ) const = 0;
};

// The implementation of this static method decides which logger
// creator will be used for the static loggers

AbstractLoggerCreator& theLoggerCreator();

} /* end namespace logging */

} /* end namespace scai */
