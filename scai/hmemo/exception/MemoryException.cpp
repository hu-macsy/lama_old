/**
 * @file MemoryException.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of methods for MemoryException
 * @author Thomas Brandes
 * @date 16.10.2015
 */

#include <scai/hmemo/exception/MemoryException.hpp>

namespace scai
{

namespace hmemo
{

MemoryException::MemoryException()
{
    mMessage = "MemoryException";
}

MemoryException::MemoryException( const std::string& message )
    : mMessage( message )
{
    mMessage += " @Memory";
}

MemoryException::~MemoryException() throw()
{
}

const char* MemoryException::what() const throw ()
{
    return mMessage.c_str();
}

} /* end hmemo */

} /* end scai */
