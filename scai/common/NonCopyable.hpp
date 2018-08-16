/**
 * @file common/NonCopyable.hpp
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
 * @brief Definition of help class that disables default copy constructors.
 * @author Jiri Kraus
 * @date 04.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

namespace scai
{

namespace common
{

/** @brief Base class to disable compiler generated copy constructor and assignment operator.
 *
 *  All classes where the default copy constructor might result in serious problems should
 *  derive from this class.
 *
 *  In the following class dynamically allocated data is freed in the destructor. If an
 *  instance variable has been copied the destructor is called twice and the data
 *  would be freed twice.
 *
 *  \code
 *  class Example : common::NonCopyable
 *  {
 *      ~Example()
 *      {
 *         if ( mData != NULL) delete mData;
 *         mData = NULL;
 *      }
 *      private :
 *         Data* mData;
 *  };
 *  \endcode
 *
 *  A class that is derived from this base class can never be used itself in a container class.
 *
 *  \code
 *  std::vector<Example> myExamples;   // not possible
 *  std::vector<*Example> myExamples;  // allocate variables always by own constructor
 *  \endcode
 *
 * The typical error message for using it in a container class is as follows:
 *
 * NonCopyable.hpp: error »const common::NonCopyable& common::NonCopyable::operator=(const common::NonCopyable&)« is private
 */

class COMMON_DLL_IMPORTEXPORT NonCopyable
{
protected:

    NonCopyable()
    {
    }

    ~NonCopyable()
    {
    }

private:

    NonCopyable( const NonCopyable& other );

    const NonCopyable& operator=( const NonCopyable& other );
};

} /* end namespace common */

} /* end namespace scai */
