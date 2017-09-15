/**
 * @file BFError.hpp
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
 * @brief BFError.h
 * @author Jiri Kraus
 * @date 06.04.2011
 */
/*
 * BFError.h
 *
 *  Created on: 01.02.2011
 *      Author: rrehrman
 */

#pragma once

#include <string>
#include <exception>

#include <scai/common/config.hpp>

namespace scai
{

namespace benchmark
{

class COMMON_DLL_IMPORTEXPORT BFError: public std::exception
{
public:
    /**
     * @brief The default constructor to create an Exception without message.
     */
    BFError();
    /**
     * @brief The constructor to create an Exception with the given message.
     * @param[in] message   The message of the Exception.
     */
    BFError( const std::string& message );
    /**
     * @brief The destructor of the Exception.
     */
    virtual ~BFError() throw ();

    /**
     * @brief Returns the message of the Exception.
     * @return The message of the Exception.
     */
    virtual const char* what() const throw ();
private:
    /** The message of the Exception. */
    std::string m_message;
};

} // namespace benchmark

} // namespace scai
