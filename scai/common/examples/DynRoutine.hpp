/**
 * @file DynRoutine.hpp
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
 * @brief Definition of a base class inclusive factory for example of using a library module.
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#include <scai/common/Factory.hpp>

/** Base class that provides by deriving from Factory a factory with a create routine.
 *
 *  The input value for the create routine is a string, the output value a pointer to
 *  a new created DynRoutine object.
 *
 *  \code
 *       *DynRoutine create( std::string )
 *  \endcode
 */

class DynRoutine  : public scai::common::Factory<std::string, DynRoutine*>
{
public:

    /** This routine must be provided by all derived classes. */

    virtual void doIt() = 0;

    virtual ~DynRoutine()
    {
    }
};
