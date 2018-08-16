/**
 * @file UnsupportedException.hpp
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

#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

namespace scai
{

namespace common
{

/** Derived exception class that is used when an unsupported feature should result in error. */

class COMMON_DLL_IMPORTEXPORT UnsupportedException : public Exception
{
public:

    /** Enumeration type for the different handlings of unsupported features. */

    enum UnsupportedType
    {
        UNSUPPORTED_WARN,      //!< give a warning if an unsupported feature is used
        UNSUPPORTED_ERROR,     //!< throw an Unspupprted exception
        UNSUPPORTED_IGNORE,    //!< neither warning nor error
        UNSUPPORTED_UNDEFINED  //!< for convenience
    };

    UnsupportedException();

    UnsupportedException( const std::string& message );

    virtual ~UnsupportedException() throw();

    virtual const char* what() const throw();

    static UnsupportedType getUnsupportedSetting();

    /** only for test purpose this routine can be used to re-read the environment variable */

    static void resetSetting();

protected:

    std::string mMessage;

private:

    static UnsupportedType unsupportedSetting;
};

} /* end namespace common */

} /* end namespace scai */
