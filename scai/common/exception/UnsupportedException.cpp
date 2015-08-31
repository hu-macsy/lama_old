/**
 * @file LAMAInterface.hpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Eric Schricker
 * @date 10.08.2015
 * @since 2.0.0
 */

#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/Settings.hpp>

namespace scai
{

namespace common
{
    
UnsupportedException::UnsupportedType UnsupportedException::unsupportedSetting = UnsupportedException::UNSUPPORTED_UNDEFINED;

UnsupportedException::UnsupportedException()
{
}

UnsupportedException::UnsupportedException( const std::string& message )
    : mMessage( message )
{
}

UnsupportedException::~UnsupportedException() throw()
{
}

const char* UnsupportedException::what() const throw ()
{
    return mMessage.c_str();
}

UnsupportedException::UnsupportedType UnsupportedException::getUnsupportedSetting()
{
    if( unsupportedSetting == UNSUPPORTED_UNDEFINED )
    {
        std::string val = "WARN";

        //bool isSet = Settings::getEnvironment( val, "SCAI_UNSUPPORTED" );
        
        Settings::getEnvironment( val, "SCAI_UNSUPPORTED" );

        // transform to uppercase

        for( std::string::iterator p = val.begin(); val.end() != p; ++p )
        {
            *p = static_cast<char>( toupper( *p ) );
        }

//        SCAI_LOG_INFO( logger, "SCAI_UNSUPPORTED=" << val << ", setting used for LAMA" )

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
            //SCAI_LOG_ERROR( logger, "SCAI_UNSUPPORTED=" << val << ", illegal value, take WARN" )
            unsupportedSetting = UNSUPPORTED_WARN;
        }
    }

    return unsupportedSetting;
}

} /* end common */

} /* end scai */
