/**
 * @file GPIException.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief GPIException.cpp
 * @author Lauretta Schubert
 * @date 25.02.2014
 */

#include <scai/dmemo/gpi/GPIException.hpp>

#include <sstream>

namespace scai
{

namespace dmemo
{

GPIException::GPIException( const std::string& message, const int gpiStatus )
{
    std::ostringstream oss;
    oss << message << " Cause: " << gpiStatus;
    mMessage = oss.str();
}

GPIException::~GPIException() throw ()
{
}

} /* end namespace dmemo */

} /* end namespace scai */

