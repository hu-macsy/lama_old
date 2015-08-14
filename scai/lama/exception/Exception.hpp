/**
 * @file Exception.hpp
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
 * @brief Interface of the class Exception.
 * @author Jiri Kraus
 * @date 02.03.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/Exception.hpp>

// logging
#include <scai/logging.hpp>

#include <exception>
#include <string>

/**
 * @brief The namespace lama holds everything of the C++ Library lama
 */
namespace lama
{

/**
 * @brief The class Exception represents a general exception in lama.
 */
class COMMON_DLL_IMPORTEXPORT Exception1: public std::exception
{
public:

    enum UnsupportedType
    {
        UNSUPPORTED_WARN, UNSUPPORTED_ERROR, UNSUPPORTED_IGNORE, UNSUPPORTED_UNDEFINED
    };

    /**
     * @brief The default constructor creates an Exception with no message.
     */
    Exception1();
    /**
     * @brief This constructor creates an Exception with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
    Exception1( const std::string& message );
    /**
     * @brief The destructor destroys this Exception.
     */
    virtual ~Exception1() throw ();
    /**
     * @brief what() returns the message of this Exception.
     *
     * @return the message of this Exception.
     */
    virtual const char* what() const throw ();

    /**
     *  This method returns the actual setting for how to deal with
     *  unsupported features.
     */
    static UnsupportedType getUnsupportedSetting();

protected:

    std::string mMessage;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private    :

    static UnsupportedType unsupportedSetting;
};

}
//namespace lama

#define LAMA_WEG_THROWEXCEPTION( msg )                                                 \
    {                                                                              \
        std::ostringstream errorStr;                                               \
        errorStr<<"Exception in line "<<__LINE__<<" of file "<<__FILE__<<"\n";     \
        errorStr<<"    Message: "<<msg<<"\n";                                      \
        common::Exception::addCallStack( errorStr );                               \
        throw common::Exception( errorStr.str() );                                 \
    }

/** This macro should be used to give hints about unsupported features.
 *
 *  It should only be used in cases where less performant solutions are still
 *  available.
 *
 *  By setting the environment variable LAMA_UNSUPPORTED to WARN only warnings
 *  will be given. For IGNORE no message is given at all. Otherwise an exception
 *  is thrown.
 */
#define LAMA_UNSUPPORTED( msg )                                                    \
    {                                                                              \
        if ( lama::Exception1::getUnsupportedSetting() !=                           \
                lama::Exception1::UNSUPPORTED_IGNORE )                              \
        {                                                                          \
            std::ostringstream errorStr;                                           \
            errorStr << "Unsupported at line ";                                    \
            errorStr << __LINE__ << " of file " << __FILE__ << "\n";               \
            errorStr << "    Message: " << msg << std::endl;                       \
            errorStr << "Use environment variable LAMA_UNSUPPORTED";               \
            errorStr << " (WARN or IGNORE) to get rid of this message";            \
            errorStr << std::endl;                                                 \
            if ( lama::Exception1::getUnsupportedSetting() ==                       \
                    lama::Exception1::UNSUPPORTED_ERROR )                           \
            {                                                                      \
                throw common::Exception( errorStr.str() );                         \
            }                                                                      \
            else                                                                   \
            {                                                                      \
                std::cout << errorStr.str() << std::endl;                          \
            }                                                                      \
        }                                                                          \
    }

