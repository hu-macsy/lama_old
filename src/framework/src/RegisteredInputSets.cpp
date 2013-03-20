/**
 * @file RegisteredInputSets.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief RegisteredInputSets.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
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

#include <framework/src/BenchmarkPrinter.hpp>
#include <framework/src/benchmark_framework.hpp>
#include <framework/src/frame_stdlib.hpp>

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
