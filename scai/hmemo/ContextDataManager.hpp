/**
 * @file ContextDataManager.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of a new dynamic array class where the array data can be
 *        used in different contexts and where the data is moved implicitly
 * @author Thomas Brandes
 * @date 19.07.2015
 */

#pragma once

// base classes
#include <scai/common/Printable.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/AccessKind.hpp>

// local libray
#include <scai/hmemo/ContextData.hpp>
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

// std
#include <vector>
#include <memory>

#define MEMORY_MAX_CONTEXTS 4

namespace scai
{

namespace hmemo
{

/* References to ContextData is needed for the lifetime of accesses.
 * ContextData& is not useful as ContextData can move in dynamic arrays.
 * Therefore the index to the array of ContextData is used.
 */

typedef size_t ContextDataIndex;

/** An object of this class manages the different incarnations of an array
 *  at different contexts.
 */

class COMMON_DLL_IMPORTEXPORT ContextDataManager : private common::NonCopyable, public common::Printable
{
public:

    ContextDataManager();

    ~ContextDataManager();

    ContextDataManager( ContextDataManager&& other ) noexcept;

    ContextDataManager& operator=( ContextDataManager&& other );

    /** Get the context data for a given context. A new entry can be created.
     *  This routine does not any locks or handling of data allocation or transfers.
     */

    ContextDataIndex getContextData( ContextPtr context );

    /** Get the data for a given memory. A new entry can be created.
     *  This routine does not any locks or handling of data allocation or transfers.
     */

    ContextDataIndex getMemoryData( MemoryPtr context );

    /** This routine provides context data for an access at a given context.
     *
     *  @param[in] context context at which data is needed
     *  @param[in] kind    kind of access, read or write
     *  @param[in] allocSize number of bytes needed for the array data in case of insufficient data
     *  @param[in] validSize number of bytes to be copied from a valid location
     *
     *  @returns  index to a corresponding entry for ContextData.
     */

    ContextDataIndex acquireAccess( ContextPtr context, common::AccessKind kind, size_t allocSize, size_t validSize );

    /** This routine must be called when an access is released, otherwise further accesses are not allowed.
     *
     *  @param[in] index   index to a corresponding entry for ContextData
     *  @param[in] kind    kind of access, read or write
     */

    void releaseAccess( ContextDataIndex index, common::AccessKind kind );

    /**
     * @brief Query the capacity ( in number of elements ) at a certain context.
     */
    size_t capacity( ContextPtr context ) const;

    /**
     * @brief Query if data is valid in a certain context
     */
    bool isValid( ContextPtr context ) const;

    /** Wait for last outstanding memory transfer. */

    void wait();

    /** Return number of all accesses. */

    int locked() const;

    /** Return number of accesses for a certain kind. */

    int locked( common::AccessKind kind ) const;

    void invalidateAll();

    /** Frees all allocated data but keeps records. */

    void purge();

    /** Copy operator */

    void copyAllValidEntries( const ContextDataManager& other, const size_t size );

    /** Operator [] gives access to the ContextData by a reference. */

    ContextData& operator[]( ContextDataIndex ref );

    ContextData& operator[]( ContextPtr context );

    ContextData& operator[]( MemoryPtr context );

    /** Swap of ContextDataManager required for swap of LAMA arrays. */

    void swap( ContextDataManager& other );

    /** prefetch: starts memory transfer to context asynchronously if valid data is required. */

    void prefetch( ContextPtr context, size_t size );

    void setValidData( ContextPtr context, const ContextDataManager& other, const size_t size );

    /** This routine tries to find a context where valid data is available */

    ContextPtr getValidContext( const ContextPtr prefContext ) const;

    ContextPtr getFirstTouchContextPtr() const;

    void reserve( ContextPtr context, const size_t size, const size_t validSize );

    void resize( const size_t size, const size_t validSize );

    void init( const void* data, const size_t size );

    virtual void writeAt( std::ostream& stream ) const;

    bool isInitialized() const
    {
        return findValidData() < mContextData.size();
    }

    /** Copy 'some' data from a valid context to the host */

    void getData( void* data, const size_t offset, const size_t size );

    /** Copy 'some' data from the host to the valid context */

    void setData( const void* data, const size_t offset, const size_t dataSize, const size_t allocSize );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Return context data that is valid currently.
     *
     *  This routine throws an exception if there no valid data at all.
     *  This might happen for zero-sized arrays (do not call it then at all)
     *  or for uninitialized arrays (should give a warning).
     */

    ContextDataIndex findValidData() const;

private:

    std::vector<ContextData> mContextData; // Incarnations of the array at different contexts

    std::unique_ptr<tasking::SyncToken> mSyncToken; //!<  outstanding transfers

    ContextDataIndex findContextData( ContextPtr context ) const;

    ContextDataIndex findMemoryData( MemoryPtr memory ) const;

    /** Help routine that copies valid data from one context to the other.
     *
     *  @param[in] source specifies the valid data that is copied
     *  @param[in] size is the number of bytes to copy
     *  @param[out] target specifies the new place to which the data is transferred.
     *
     *  The source and target data can belong to different arrays so this routine
     *  is also used for copy or assignment operations.
     *
     *  If a direct transfer from source to target is not possible (unsupported by
     *  the context) it will be tried to copy the data by involving a tempoarary
     *  copy in the Host context (copy belongs to the array with this ContextDataManager).
     */

    void fetch( ContextData& target, const ContextData& source, size_t size );

    tasking::SyncToken* fetchAsync( ContextData& target, const ContextData& source, size_t size );

    mutable std::mutex mAccessMutex; // needed to make accesses thread-safe, must not be recursive
    mutable std::condition_variable_any mAccessCondition;  // notify if all accesses are released

    void lockAccess( common::AccessKind kind, ContextPtr context );

    void unlockAccess( common::AccessKind kind );

    int mLock[static_cast<int>( common::AccessKind::MaxAccessKind )];

    ContextPtr accessContext;    // context that has locked

    std::thread::id accessThread;     // thread that has locked

    bool  multiContext;   // multiple reads at different Context
    bool  multiThreaded;  // multiple reads by different threads

    bool hasAccessConflict( common::AccessKind kind ) const;
};


} /* end namespace hmemo */

} /* end namespace scai */
