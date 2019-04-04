/**
 * @file ContextData.hpp
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
 * @brief Records for allocated data in a certain context.
 * @author Thomas Brandes
 * @date 14.07.2011
 */
#pragma once

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/hmemo/Memory.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/config.hpp>
#include <scai/common/macros/throw.hpp>

namespace scai
{

namespace tasking
{

class SyncToken;

}

namespace hmemo
{

/** @brief Objects of this class are used to manage different incarnations
 *         of a LAMA array at different contexts.
 *
 *  As this class might be used in a container, default constructor
 *  and copy constructor are provided.
 *
 *  The destructor will not free the allocated data as this might result
 *  in memory corruptions. The free method must be called explicitly to
 *  free the allocated data. This might be done by other constructors.
 */

class COMMON_DLL_IMPORTEXPORT ContextData: public common::Printable
{
private:

    size_t size; //!<  allocated size stands also for capacity

    MemoryPtr mMemory; //!<  shared pointer to the context

    void* pointer; //!<  pointer to the data on the context

    bool valid;     //!<  is true if data at context is valid

    bool allocated; //!<  is true if data has been allocated by context

public:

    void* get()
    {
        return pointer;
    }

    MemoryPtr getMemoryPtr() const
    {
        return mMemory;
    }

    const Memory& getMemory() const
    {
        return *mMemory;
    }

    const void* get() const
    {
        return pointer;
    }

    /** Constructor, memory object must always be given. */

    ContextData( MemoryPtr memory );

    /** Move constructor, allows use of this class container classes.
     *  The clause noexcept is mandatory and also safe here.
     */
    ContextData( ContextData&& other ) noexcept;

    /** Destructor, will also call free if not done before */

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

    /** Make reservation of enough memory.
     *
     *  @param newSize    number of entries for new allocated memory
     *  @param validSize  number of values that must be copied
     *  @param inUse      throw exception if reallocate is required
     */
    void reserve( const size_t newSize, const size_t validSize, bool inUse );

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

    tasking::SyncToken* copyFromAsync( const ContextData& source, size_t size );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    // disable the default copy constructor

    ContextData( const ContextData& ) = delete;

    ContextData& operator=( ContextData&& other ) noexcept;
};

} /* end namespace hmemo */

} /* end namespace scai */
