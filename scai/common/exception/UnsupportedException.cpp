/**
 * @file UnsupportedException.cpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Eric Schricker
 * @date 10.08.2015
 */

// hpp
#include <scai/common/exception/UnsupportedException.hpp>

// local library
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
    if ( unsupportedSetting == UNSUPPORTED_UNDEFINED )
    {
        std::string val = "WARN";
        //bool isSet = Settings::getEnvironment( val, "SCAI_UNSUPPORTED" );
        Settings::getEnvironment( val, "SCAI_UNSUPPORTED" );

        // transform to uppercase

        for ( std::string::iterator p = val.begin(); val.end() != p; ++p )
        {
            *p = static_cast<char>( toupper( *p ) );
        }

//        SCAI_LOG_INFO( logger, "SCAI_UNSUPPORTED=" << val << ", setting used for LAMA" )

        if ( "IGNORE" == val )
        {
            unsupportedSetting = UNSUPPORTED_IGNORE;
        }
        else if ( "WARN" == val )
        {
            unsupportedSetting = UNSUPPORTED_WARN;
        }
        else if ( "ERROR" == val )
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

void UnsupportedException::resetSetting()
{
    unsupportedSetting = UNSUPPORTED_UNDEFINED;
}

} /* end common */

} /* end scai */
