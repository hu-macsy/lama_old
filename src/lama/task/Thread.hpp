/**
 * @file Thread.hpp
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
 * @brief Definition of the class Thread
 * @author Jiri Kraus
 * @date 30.03.2012
 * @since 1.0.0
 */
#ifndef LAMA_THREAD_HPP_
#define LAMA_THREAD_HPP_

// for dll_import
#include <common/config.hpp>

// boost
#include <boost/thread.hpp>

namespace lama
{

#ifdef WIN32
#define LAMA_USE_BOOST_THREADID
#endif

#ifndef LAMA_USE_BOOST_THREADID
/** boost::thread::id is only available with Boost 1.35 and higher */
#include <boost/version.hpp>

#define LAMA_BOOST_VERSION_PROVIDES_ID 103501

#if BOOST_VERSION < LAMA_BOOST_VERSION_PROVIDES_ID
#include <pthread.h>
#else
#define LAMA_USE_BOOST_THREADID
#endif
#endif

class COMMON_DLL_IMPORTEXPORT Thread
{
public:

#ifdef LAMA_USE_BOOST_THREADID
    typedef boost::thread::id Id;
#else
    typedef pthread_t Id;
#endif

    /** returns the id of the calling Thread. */
    static Id getSelf();
};

} // namespace lama

#endif // LAMA_THREAD_HPP_
