/**
 * @file common/NonCopyable.hpp
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
 * @brief Definition of help class that disables default copy constructors.
 * @author Jiri Kraus
 * @date 04.04.2011
 * @since 1.0.0
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
