/**
 * @file RegisteredInputSets.cpp
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
 * @brief RegisteredInputSets.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 */
/*
 * RegisteredInputSets.cpp
 *
 *  Created on: 31.01.2011
 *      Author: rrehrman
 */

#include <string>
#include <vector>
#include <map>

#include <scai/benchmark.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>
#include <scai/benchmark/frame_stdlib.hpp>

using namespace scai;

extern "C" bf::BaseInputSetRegistry* getInputSetRegistry();

int main( void )
{
    std::vector<std::string> files;

    try
    {
        bf::getSharedLibraries( files );
    }
    catch( std::exception& e )
    {
        bf::BenchmarkPrinter::error( e.what() );
        return 1;
    }

    bf::BaseInputSetRegistry* registry = NULL;
    std::map<std::string,std::string> inputSets;

    // try dlopen on the shared libraries, saving all registered InputSets from
    // the shared libraries in a vector of strings.
    for( unsigned int i = 0; i < files.size(); ++i )
    {
        typedef bf::BaseInputSetRegistry* (*registry_t)();
        registry_t reg_handle = NULL;

        LAMA_LIB_HANDLE_TYPE handle;

        int error = bf::loadLibAndGetFunctionHandle( reg_handle, handle, files[i].c_str(), "getInputSetRegistry" );

        if( error != 0 )
        {
            continue;
        }

        typedef void (*free_registry_t)();
        free_registry_t reg_free_handle;

        error = bf::getFunctionHandle( reg_free_handle, handle, "releaseInputSetRegistry" );
        if( error != 0 )
        {
            std::stringstream message;
            message << files[i] << " does not define releaseInputSetRegistry skipping it.";
            bf::BenchmarkPrinter::warning( message.str() );
            //bf::freeLibHandle( handle );
            continue;
        }

        // getting registry from library.
        registry = reg_handle();

        std::map<std::string,std::string> iSets;

        registry->getInputSetMap( iSets );
        inputSets.insert( iSets.begin(), iSets.end() );

        reg_free_handle();
        //bf::freeLibHandle( handle );
    }

    std::stringstream message;

    for( std::map<std::string,std::string>::const_iterator it = inputSets.begin(); it != inputSets.end(); ++it )
    {
        message << it->first << "%,";
    }
    std::string messagestr = message.str();

    // remove the last '%,'
    bf::BenchmarkPrinter::print( messagestr.substr( 0, messagestr.size() - 2 ) );
}
