/**
 * @file InvalidArgumentException.hpp
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
 * @brief Definition of exception class that is thrown when an illegal argument occurs
 * @author Thomas Brandes
 * @date 15.01.2018
 */
#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

// std
#include <string>

namespace scai
{

namespace common
{

/** Derived exception class for an exception that is thrown if an invalid argument is used. */

class COMMON_DLL_IMPORTEXPORT InvalidArgumentException : public Exception
{
public:

    /**
     * @brief This constructor creates an InvalidArgumentException with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
    inline InvalidArgumentException( const std::string& message );

    /**
     * @brief The destructor destroys this InvalidArgumentException.
     */
    inline virtual ~InvalidArgumentException() throw ();
};

InvalidArgumentException::InvalidArgumentException( const std::string& message ) :

    Exception( message )
{
    mMessage += "@InvalidArgumentException";
}

InvalidArgumentException::~InvalidArgumentException() throw ()
{
}

} /* end namespace common */

} /* end namespace scai */

