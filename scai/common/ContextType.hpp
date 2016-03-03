/**
 * @file ContextType.hpp
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
 * @brief Definition of enum type for scalar value types typically used in numerical applications.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace common
{

/** Own namespace for the enumeration type ContextType and its values. 
 *
 *  Note: This enumeration type was originally defined in the class Context so
 *        it did not require an own namespace. But this is no longer possible as
 *        Context is now a factory where ContextType is the input type to create
 *        a new Context.
 */

namespace context
{
    /** Enumeration type for the supported contexts. The type is used to select
     *  the appropriate code that will be used for the computations in the context.
     *
     *  The same context type does not imply that two different contexts can use
     *  the same data. Two CUDA contexts might allocate their own data where data
     *  must be transfered explicitly.
     */
    typedef enum 
    {
        Host,          //!< context for cpu + main memory
        CUDA,          //!< CUDA GPU device
        OpenCL,        //!< OpenCL GPU device, currently not supported
        MIC,           //!< Intel MIC
        UserContext,   //!< can be used for a new derived Context class
        MaxContext     //!< used for dimension of ContextType arrays
    } ContextType;

    /** Enumeration type for access kind, may be read or write */

    typedef enum
    {
        Read, //!<  read access to the array, can be multiple
        Write, //!<  write access to the array, only one at a time
        MaxAccessKind //!<  internal use for dimension of arrays
    } AccessKind;

    COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ContextType& type );

    COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const AccessKind& kind );

} /* end namespace context */

} /* end namespace common */

} /* end namespace scai */
