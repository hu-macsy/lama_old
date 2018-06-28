/**
 * @file AssertException.cpp
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
 * @brief Implementaion of methods for class AssertException.
 * @author Eric Schricker
 * @date 31.08.2015
 */

#include <scai/common/exception/AssertException.hpp>

namespace scai
{

namespace common
{

AssertException::AssertException() :
    mMessage( "AssertException" )
{
}

AssertException::AssertException( const std::string& message ) :
    mMessage( message )
{
    mMessage += "@AssertException";
}

AssertException::~AssertException() throw ()
{
}

const char* AssertException::what() const throw ()
{
    return mMessage.c_str();
}

} /* end namespace common */

} /* end namespace scai */


