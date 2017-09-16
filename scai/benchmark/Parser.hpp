/**
 * @file Parser.hpp
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
 * @brief Parsing of strings to create dynamic objects
 * @author Thomas Brandes
 * @date 19.09.2017
 */

#pragma once

#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>

#include <string>
#include <iostream>


namespace 
{

static inline void tokenize ( std::vector<std::string>& out, const std::string& in, const std::string& delimiters ) 
{
    scai::common::Settings::tokenize( out, in, delimiters );
}

/** @brief parse string like "name ( args )" to components name and args without blanks
 *
 *  This is used to translate BenchmarkId( benchmarkArgs ) or InputSetId( inputSetArgs )
 *  into corresponding create calls of the factory.
 */
static inline void parseCommand( std::string& name, std::string& args, const std::string& command )
{
    std::vector<std::string> tokens;

    tokenize( tokens, command, " ,():" );

    SCAI_ASSERT_GE_ERROR( tokens.size(), 1, "No command: " << command );

    name = tokens[0];

    if ( tokens.size() <= 1 )
    {
        args = "";
    }
    else
    {
        args = tokens[1];
        for ( size_t i = 2; i < tokens.size(); ++i )
        {
            args += ", " + tokens[i];
        }
    }

    std::cout << "#tokens = " << tokens.size()
              << ", name=<" << name << ">, args=<" << args << ">" << std::endl;

}

} 
