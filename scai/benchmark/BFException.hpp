/**
 * @file BFException.hpp
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
 * @brief BFException.h
 * @author jiri
 * @date 06.04.2011
 */
/**
 * @file BFException.h
 * @author jiri
 * Created on: 11.05.2010
 */
#pragma once

#include <scai/common/config.hpp>
#include <scai/benchmark/BFError.hpp>

namespace scai
{

namespace bf
{

class COMMON_DLL_IMPORTEXPORT BFException: public BFError
{
public:
    /**
     * @brief Default constructor.
     */
    BFException();
    /**
     * @brief Constructor initializes BFException with given message.
     * @param[in] message The message of the Error.
     */
    BFException( const std::string& message );
    /**
     * @brief Destructor.
     */
    virtual ~BFException() throw ();
};

}

}

#ifndef LAMA_CHECK_BENCHMARK
/**
 * Checks an expression and prints out a warning at compile time, if
 * expression was not true.
 */
#define LAMA_CHECK_BENCHMARK(exp)                                              \
    if (!(exp))                                                                     \
    {                                                                               \
        std::ostringstream errorStr;                                                \
        errorStr<<"Warning: Incorrect Results.\n";                                  \
        throw bf::BFException(errorStr.str( ));                                \
    }

#endif
