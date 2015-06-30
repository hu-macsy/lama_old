/**
 * @file Context.hpp
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
 * @brief Records for allocated data in a certain context.
 *
 * @author Thomas Brandes
 * @date 14.07.2011
 */
#pragma once

namespace memory
{

struct COMMON_DLL_IMPORTEXPORT ContextData: private common::NonCopyable
{
    enum AccessKind
    {
        Read, //!<  read access to the array, can be multiple
        Write, //!<  write access to the array, only one at a time
        MaxAccessKind //!<  internal use for dimension of arrays
    };

    ContextPtr context; //!<  shared pointer to the context

    void* pointer; //!<  pointer to the data on the context

    size_t size; //!<  size of a single element

    bool allocated; //!<  is true if data has been allocated by context

    bool valid; //!<  is true if there is a valid copy on the context

    unsigned char lock[MaxAccessKind]; //!<  read, write lock

    /** Constructor, context must always be given. */

    ContextData( ContextPtr context );

    ~ContextData();

    /** allocate data for the Context array on the context */

    void allocate( const size_t size );

    /** set reference instead of own data allocation. */

    void setRef( void* reference, const size_t size );

    /** Reallocate new memory on the context
     *
     *  @param newSize   number of entries for new allocated memory
     *  @param saveSize  number of entries to copy back from old memory
     *
     */
    void realloc( const size_t newSize, const size_t saveSize );

    /** free data for the Context array on the context */

    void free();

    bool isPinned() const;

    void setPinned() const;

    void setCleanFunction( boost::function<void( void* )> cleanFunktion ) const;

private:

    mutable bool pinned;

    mutable boost::function<void( void* )> mCleanFunktion;

    ContextData(); // disable default constructor
};

}

