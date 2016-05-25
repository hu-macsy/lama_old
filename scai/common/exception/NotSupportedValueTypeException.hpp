/**
 * @file NotSupportedValueTypeException.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Definition of derived exception class for unsupported values.
 * @author Eric Schricker
 * @date 12.11.2015
 */

// base class
#include <scai/common/exception/Exception.hpp>

// std
#include <string>

namespace scai
{

namespace common
{

/** Derived exception class that is used to throw an exception in all cases
 *  where an illegal or unsupported value is used. 
 */

class COMMON_DLL_IMPORTEXPORT NotSupportedValueTypeException : public Exception
{
public:
    /**
     * @brief The default constructor creates an NotSupportedValueTypeException with no message.
     */
	NotSupportedValueTypeException( );

    /**
     * @brief This constructor creates an NotSupportedValueTypeException with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
	NotSupportedValueTypeException( const std::string& message );

    /**
     * @brief The destructor destroys this NotSupportedValueTypeException.
     */
    virtual ~NotSupportedValueTypeException() throw ();

    virtual const char* what() const throw();

protected:

    std::string mMessage;
};

} /* end namespace common */

} /* end namespace scai */

