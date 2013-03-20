/**
 * @file LAMAArguments.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief LAMAArguments.hpp
 * @author Michael Drost
 * @date 07.03.2012
 * $Id$
 */
#ifndef LAMA_LAMAINIT_HPP_
#define LAMA_LAMAINIT_HPP_

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <LAMAInstance.hpp>

namespace lama
{

class LAMAArguments
{
    friend class LAMAInstance;
    int mArgc;
    char** mArgv;
//    static std::vector<LAMAArgumentFilter> mFilters;
    LAMAArguments( int argc, char** argv );
    static std::auto_ptr<LAMAArguments> sInstance;
public:
    static void init( int argc, char** argv, bool keepArguments = false );
    static char** getArgv();
    static int getArgc();
    static void finalize();
    static bool isInitialized();
    virtual ~LAMAArguments();
};

}
#endif // LAMA_LAMAINIT_HPP_
