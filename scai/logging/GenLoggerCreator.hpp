/**
 * @file GenLoggerCreator.hpp
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
 * @brief Definition of a specific Logger creator that implements the abstract one.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// base classes
#include <scai/logging/AbstractLoggerCreator.hpp>

namespace scai
{

namespace logging
{

class GenLoggerCreator: public AbstractLoggerCreator
{

public:

    virtual ~GenLoggerCreator();

    /**
     * @brief Function to create a new instance at a given level.
     * @param[in] name identifies the instance must not contain any dots
     * @param parent is existent parent logger
     * @return Pointer to the new created instance
     *
     * Note: the new created instance will be added as a new son to parent.
     */

    virtual Logger* create( const std::string& name, Logger* parent ) const;

    /** Getter for the root instance of all loggers. */

    virtual Logger& getRoot() const;

    /** Get the single instantiation of this class */

    static GenLoggerCreator& getTheCreator();

private:

    /** Only one instance of this object is generated.
     *
     */
    GenLoggerCreator()
    {
    }

    static GenLoggerCreator* theCreator; //!< singleton instantiation
};

} /* end namespace logging */

} /* end namespace scai */
