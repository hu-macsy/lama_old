/**
 * @file CUDAException.hpp
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
 * @brief Definition of exception class that is thrown when an CUDA problem occurs.
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

/** Derived exception class for an exception that is thrown if any CUDA problem occurs. */

class COMMON_DLL_IMPORTEXPORT CUDAException : public Exception
{
public:

    /**
     * @brief This constructor creates an CUDAException with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
    inline CUDAException( const std::string& message );

    /**
     * @brief The destructor destroys this CUDAException.
     */
    inline virtual ~CUDAException() throw ();
};

CUDAException::CUDAException( const std::string& message ) :

    Exception( message )
{
    mMessage += "@CUDAException";
}

CUDAException::~CUDAException() throw ()
{
}

} /* end namespace common */

} /* end namespace scai */

