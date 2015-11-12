/**
 * @file MemoryException.hpp
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
 * @brief Derived exception class for exceptions in memory allocation/deallocation
 * @author Thomas Brandes
 * @date 10.10.2015
 */

#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

namespace scai
{

namespace hmemo
{

/** Derived class needed to catch exception for memory problems */

class COMMON_DLL_IMPORTEXPORT MemoryException : public scai::common::Exception
{
public:

    MemoryException();
    
    MemoryException( const std::string& message );

    virtual ~MemoryException() throw();
    
    virtual const char* what() const throw();
    
protected:
        
    std::string mMessage;
};

} /* end namespace common */

} /* end namespace scai */
