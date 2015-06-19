/**
 * @file lama/NonCopyable.hpp
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
#ifndef LAMA_NONCOPYABLE_HPP_
#define LAMA_NONCOPYABLE_HPP_

// for dll_import
#include <common/config.hpp>

namespace lama
{

/** Base class to disable compiler generated copy constructor and assignment operator. */
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

}

#endif // LAMA_NONCOPYABLE_HPP_
