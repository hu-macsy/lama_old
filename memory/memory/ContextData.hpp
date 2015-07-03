/**
 * @file ContextData.hpp
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

#include <logging/logging.hpp>

#include <common/config.hpp>
#include <common/Exception.hpp>
#include <common/Printable.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

namespace memory
{

class Context;
class SyncToken;

/** Context pointers will be always const, so context can never be modified. */

typedef boost::shared_ptr<const Context> ContextPtr;

/** @brief Objects of this class are used to manage different incarnations
 *         of a LAMA array at different contexts.
 *
 *  As this class might be used in a container, default constructor
 *  and copy constructor are provided.
 *
 *  The destructor will not free the allocated data as this might result
 *  in memory corruptions. The free method must be called explicitly to
 *  free the allocated data.
 */

class COMMON_DLL_IMPORTEXPORT ContextData: public Printable
{
private:

    size_t size; //!<  allocated size stands also for capacity

    ContextPtr mContext; //!<  shared pointer to the context

    void* pointer; //!<  pointer to the data on the context

    bool valid;     //!<  is true if data at context is valid

    bool allocated; //!<  is true if data has been allocated by context

public:

    enum AccessKind
    {
        Read, //!<  read access to the array, can be multiple
        Write, //!<  write access to the array, only one at a time
        MaxAccessKind //!<  internal use for dimension of arrays
    };

    void* get() { return pointer; }

    ContextPtr context() const { return mContext; }

    const void* get() const { return pointer; }

    /** Constructor, context must always be given. */

    ContextData( ContextPtr context );

    ContextData();  // allow default constructor for container

    /** Destructor, will NOT free allocated data at a context. */

    ~ContextData();

    size_t capacity() const 
    {
        return size;
    }

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

    void reserve( const size_t newSize, const size_t validSize );

    /** free allocated data at the corresponding context */

    void free();

    void setValid( bool flag ) 
    {
        valid = flag;
    }

    bool isValid() const
    {
        return valid;
    }

    virtual void writeAt( std::ostream& stream ) const;

    void copyFrom( const ContextData& source, size_t size );

    SyncToken* copyFromAsync( const ContextData& source, size_t size );

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

};

}

