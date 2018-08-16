/**
 * @file GenLoggerCreator.cpp
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
 * @brief Implemenation of methods for class GenLoggerCreator.
 * @author Thomas Brandes
 * @date 01.03.2011
 */

// local library
#include <scai/logging.hpp>
#include <scai/logging/GenLoggerCreator.hpp>
#include <scai/logging/GenLogger.hpp>

// std
#include <iostream>
#include <cstdlib>         // import getenv
#include <cstdio>          // FILE
#include <stdexcept>       // runtime_error

#undef DEBUGGING

using namespace std;

namespace scai
{

namespace logging
{

// GenLoggerCreator becomes the 'static' global logger creator.

AbstractLoggerCreator& theLoggerCreator()
{
    return GenLoggerCreator::getTheCreator();
}

GenLoggerCreator* GenLoggerCreator::theCreator = NULL;

GenLoggerCreator& GenLoggerCreator::getTheCreator()
{
    if ( !theCreator )
    {
        theCreator = new GenLoggerCreator();
    }

    return *theCreator;
}

GenLoggerCreator::~GenLoggerCreator()
{
}

/********************************************************************
 *  GenLoggerCreator:: getRoot()                                     *
 ********************************************************************/

Logger& GenLoggerCreator::getRoot() const
{
    return GenLogger::getRoot();
}

Logger* GenLoggerCreator::create( const std::string& name, Logger* parent ) const
{
    return new GenLogger( name, parent );
}

} /* end namespace logging */

} /* end namespace scai */
