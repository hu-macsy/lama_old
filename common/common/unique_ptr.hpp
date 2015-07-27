/**
 * @file common/unique_ptr.hpp
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
 * @brief Embedding unique_ptr in common namespace.
 * @author Thomas Brandes
 * @date 15.07.2015
 */

#pragma once

#if __cplusplus > 199711L
    #include <memory>
#else
 
    // Boost 1.57 and higher provides same unique_ptr as C++
    // but unfortunately not the previous versions
    // #include <boost/interprocess/smart_ptr/unique_ptr.hpp>

    #include <boost/scoped_array.hpp>
    #include <memory>

#endif

namespace common
{
#if __cplusplus > 199711L

    using std::unique_ptr;

    /** C++11: here we can use unique_ptr for a scoped array */

    template<typename T>
    class scoped_array public std::unique_ptr<T[]>
    {
    public:
        scoped_array( T* ptr ) : std::unique_ptr<T[]>( ptr )
        {
        }
    };

#else

    using boost::scoped_array;

    template<typename T>
    class unique_ptr : public std::auto_ptr<T>
    {
    public:
        unique_ptr( T* ptr ) : std::auto_ptr<T>( ptr ) 
        {
        }
        unique_ptr() : std::auto_ptr<T>() 
        {
        }
    };
#endif
}
