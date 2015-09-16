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

#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

namespace scai
{

namespace common
{

class COMMON_DLL_IMPORTEXPORT UnsupportedException : public Exception
{
public:

    enum UnsupportedType
    {
        UNSUPPORTED_WARN, UNSUPPORTED_ERROR, UNSUPPORTED_IGNORE, UNSUPPORTED_UNDEFINED
    };
    
    UnsupportedException();
    
    UnsupportedException( const std::string& message );
    virtual ~UnsupportedException() throw();
    
    virtual const char* what() const throw();
    
    static UnsupportedType getUnsupportedSetting();
    
protected:
        
    std::string mMessage;

private:

    static UnsupportedType unsupportedSetting;
};

} /* end namespace common */

} /* end namespace scai */

/** This macro should be used to give hints about unsupported features.
 *
 *  It should only be used in cases where less performant solutions are still
 *  available.
 *
 *  By setting the environment variable COMMON_UNSUPPORTED to WARN only warnings
 *  will be given. For IGNORE no message is given at all. Otherwise an exception
 *  is thrown.
 */
 
#define SCAI_UNSUPPORTED( msg )                                                \
{                                                                              \
    if ( scai::common::UnsupportedException::getUnsupportedSetting() !=        \
            scai::common::UnsupportedException::UNSUPPORTED_IGNORE )           \
    {                                                                          \
        std::ostringstream errorStr;                                           \
        errorStr << "Unsupported at line ";                                    \
        errorStr << __LINE__ << " of file " << __FILE__ << "\n";               \
        errorStr << "    Message: " << msg << std::endl;                       \
        errorStr << "Use environment variable SCAI_UNSUPPORTED";               \
        errorStr << " (WARN or IGNORE) to get rid of this message";            \
        errorStr << std::endl;                                                 \
        if ( scai::common::UnsupportedException::getUnsupportedSetting() ==    \
        		scai::common::UnsupportedException::UNSUPPORTED_ERROR )        \
        {                                                                      \
            throw scai::common::UnsupportedException( errorStr.str() );        \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            std::cout << errorStr.str() << std::endl;                          \
        }                                                                      \
    }                                                                          \
}
