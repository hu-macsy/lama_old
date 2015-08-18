/**
 * @file ContextDataManager.hpp
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
 * @brief Definition of a new dynamic array class where the array data can be
 *        used in different contexts and where the data is moved implicitly
 *        when corresponding read/write accesses are required
 * @author Thomas Brandes
 * @date 19.07.2015
 */

#pragma once

#include <scai/hmemo/ContextData.hpp>
#include <scai/tasking/SyncToken.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/logging.hpp>

#include <scai/common/NonCopyable.hpp>

#include <vector>

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

class ContextDataManager : private common::NonCopyable, public Printable
{
public:

    ContextDataManager();

    ~ContextDataManager();

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

    ContextDataIndex acquireAccess( ContextPtr context, AccessKind kind, size_t allocSize, size_t validSize );

    /** This routine must be called when an access is released, otherwise further accesses are not allowed. 
     *
     *  @param[in] index   index to a corresponding entry for ContextData
     *  @param[in] kind    kind of access, read or write
     */

    void releaseAccess( ContextDataIndex index, AccessKind kind );

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
  
    /** Return true if there is at least one access to any context data. 
     *  For a locked array further write access is not possible.
     */

    bool locked() const;

    /** Return true if there is at least one write access to any context data. */

    bool locked( AccessKind kind ) const;

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

    ContextPtr getValidContext( const ContextType preferredType );

    void reserve( ContextPtr context, const size_t size, const size_t validSize );

    void resize( const size_t size, const size_t validSize );

    virtual void writeAt( std::ostream& stream ) const;

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

    std::auto_ptr<tasking::SyncToken> mSyncToken; //!<  outstanding transfers

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

    mutable common::Thread::Mutex mAccessMutex; // needed to make accesses thread-safe, must not be recursive 
    mutable common::Thread::Condition mAccessCondition;  // notify if all accesses are released

    void lockAccess( AccessKind kind, ContextPtr context );

    void unlockAccess( AccessKind kind );

    int mLock[context::MaxAccessKind];

    ContextPtr accessContext;    // context that has locked
    common::Thread::Id accessThread;     // thread that has locked

    bool  multiContext;   // multiple reads at different Context
    bool  multiThreaded;  // multiple reads by different threads

    bool hasAccessConflict( AccessKind kind ) const;
};


} /* end namespace hmemo */

} /* end namespace scai */
