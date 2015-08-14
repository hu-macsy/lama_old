/**
 * @file Exception.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Exception.cpp
 * @author Jiri Kraus
 * @date 01.03.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/exception/Exception.hpp>
#include <scai/lama/Settings.hpp>

#include <cstdio>
#include <sstream>

#ifndef _WIN32
#include <execinfo.h>
#endif //_WIND32

#ifdef __GNUC__
#include <cxxabi.h>
#include <cstring>
#endif // __GNUC__

namespace lama
{

SCAI_LOG_DEF_LOGGER( Exception1::logger, "Exception" )

Exception1::UnsupportedType Exception1::unsupportedSetting = Exception1::UNSUPPORTED_UNDEFINED;

Exception1::UnsupportedType Exception1::getUnsupportedSetting()
{
    if( unsupportedSetting == UNSUPPORTED_UNDEFINED )
    {
        std::string val = "WARN";

        bool isSet = Settings::getEnvironment( val, "LAMA_UNSUPPORTED" );

        if( !isSet )
        {
            SCAI_LOG_WARN( logger, "LAMA_UNSUPPORTED not set, default is WARN" )
        }

        // transform to uppercase

        for( std::string::iterator p = val.begin(); val.end() != p; ++p )
        {
            *p = static_cast<char>( toupper( *p ) );
        }

        SCAI_LOG_INFO( logger, "LAMA_UNSUPPORTED=" << val << ", setting used for LAMA" )

        if( "IGNORE" == val )
        {
            unsupportedSetting = UNSUPPORTED_IGNORE;
        }
        else if( "WARN" == val )
        {
            unsupportedSetting = UNSUPPORTED_WARN;
        }
        else if( "ERROR" == val )
        {
            unsupportedSetting = UNSUPPORTED_ERROR;
        }
        else
        {
            SCAI_LOG_ERROR( logger, "LAMA_UNSUPPORTED=" << val << ", illegal value, take WARN" )
            unsupportedSetting = UNSUPPORTED_WARN;
        }
    }

    return unsupportedSetting;
}

Exception1::Exception1()
{
}

Exception1::Exception1( const std::string& message )
    : mMessage( message )
{
    SCAI_LOG_WARN( logger, "EXCEPTION: " << message )
}

Exception1::~Exception1() throw ()
{
}

const char* Exception1::what() const throw ()
{
    return mMessage.c_str();
}

}
//namespace lama
