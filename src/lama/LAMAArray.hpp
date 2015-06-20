/**
 * @file LAMAArray.hpp
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
 * @date 14.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_LAMA_ARRAY_HPP_
#define LAMA_LAMA_ARRAY_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/Context.hpp>
#include <lama/SyncToken.hpp>
#include <lama/Scalar.hpp>

// boost
#include <boost/thread/recursive_mutex.hpp>

#include <vector>

/** Number of contexts that might be used in maximum. This number
 *  is used for reservation of entries but does not imply any restrictions.
 */

#define LAMA_MAX_CONTEXTS 4

namespace lama
{

// Forward declaration of friend classes.

template<typename ValueType>
class ReadAccess;

template<typename ValueType>
class WriteAccess;

/** Common base class for typed LAMAArray. */

class COMMON_DLL_IMPORTEXPORT _LAMAArray: public Printable
{

public:

    virtual ~_LAMAArray()
    {
    }

    /**
     * @brief Query the value type of the matrix elements, e.g. DOUBLE or FLOAT.
     */
    virtual Scalar::ScalarType getValueType() const = 0;

    /**
     * @brief Create an empty array of a certain type.
     *
     * @param type specifies the type of the array to be created.
     *
     */

    static _LAMAArray* create( const Scalar::ScalarType type );

    /**
     * @brief Query the current size of the LAMA array, i.e. number of entries.
     *
     * @return the number of entries of the array.
     */
    inline IndexType size() const;

    virtual ContextPtr getValidContext( const Context::ContextType preferredType = Context::Host ) const = 0;

protected:

    explicit _LAMAArray( const IndexType n )
                    : mSize( n ), constFlag( false )
    {
    }

    IndexType mSize; //!< number of entries for the context array, common for all contexts

    bool constFlag; //!< if true the array cannot be written

    static size_t nContextIndex; // stands for no valid index
};

/**
 * @brief LAMAArray is the base array container for all compute relevant data within LAMA.
 *
 * @tparam ValueType is the type stored in this container.
 *
 * LAMAArray its contents on all supported locations, e.g. Host, CUDA and OpenCL. It transparently handles
 * synchronization between the locations. To enforce the consistency of the data a LAMAArray can be only
 * indirectly accessed via a ReadAccess or a WriteAccess.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArray: public _LAMAArray
{
    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;

public:

    /**
     * @brief LAMAArray() creates an empty LAMAArray with size 0
     */
    LAMAArray();

    /**
     * @brief LAMAArray( const IndexType n ) creates a LAMAArray of size n
     *
     * @param[in] n the size of the LAMAArray to create
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n and allocates uninitialized Host memory.
     */
    explicit LAMAArray( const IndexType n );

    /**
     * @brief Creates a LAMAArray of size n.
     *
     * @param[in] n     the size of the LAMAArray to create
     * @param[in] value the value to initialize the container contens with
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n, allocates Host memory and fills the Host memory with
     * the passed value.
     */
    LAMAArray( const IndexType n, const ValueType& value );

    /**
     * @brief Creates a LAMAArray of size n.
     *
     * @param[in] n         the size of the LAMAArray to create
     * @param[in] values    the values to initialize the container contens with
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n, allocates Host memory and fills the Host memory with
     * the passed values.
     */
    template<typename OtherValueType>
    LAMAArray( const IndexType n, const OtherValueType* const values );

    /**
     * @brief Creates a copy of the passed LAMAArray.
     *
     * @param[in] other the LAMAArray to copy
     *
     * LAMAArray(const LAMAArray<ValueType>& other) copies the passed LAMAArray. The container contens is copied for all currently valid
     * Locations.
     */
    LAMAArray( const LAMAArray<ValueType>& other );

    /**
     * @brief Destructor, releases all used resources.
     */
    virtual ~LAMAArray();

    /**
     * @brief Copies the passed LAMAArray into this.
     *
     * @param[in] other the LAMAArray to copy
     *
     * The assignment operator copies the passed LAMAArray.
     * The container content is copied for all contexts where a
     * valid copy is available (at least one).
     */
    LAMAArray<ValueType>& operator=( const LAMAArray<ValueType>& other );

    /*
     * @brief Checks if the LAMAArray referenced by other is the same as the LAMAArray reference by this.
     *
     * @param[in]   other   the LAMAArrayC to compare this with.
     * @return              if this and other are referencing the same LAMAArray.
     */
    bool operator==( const LAMAArray<ValueType>& other ) const
    {
        return &other == this;
    }

    bool operator!=( const LAMAArray<ValueType>& other ) const
    {
        return &other != this;
    }

    /**
     * @brief Copies the passed LAMAArray into this.
     *
     * @param[in] other     the LAMAArray to copy
     * @param[in] context   the context where the assignment should be carried out
     *
     * The assign method copies the passed LAMAArray.
     * The container content is copied for the passed contexts. If necessary other is
     * copied to context to carry this out.
     */
    void assign( const LAMAArray<ValueType>& other, ContextPtr context );

    /**
     * @brief Swaps the contens of this with other.
     *
     * @param[in] other the LAMAArray to swap the contens with.
     */
    void swap( LAMAArray<ValueType>& other );

    /**
     * @brief Prefetches the contents of the container to the passed context.
     *
     * @param[in] context  the context to prefetch to
     *
     * This method prefetches the contents of the container to the context.
     * If this valid at location nothing happens,if not a transfer from a valid location
     * to the passed location is started. Because the transfer is handled by LAMAArray and to
     * maintain the consistency of the container only one running transfer can exist at any
     * point in time. There for if two prefetches to two different invalid locations are
     * started one after the other the second transfer does not start before the first one is finished.
     */
    void prefetch( ContextPtr context ) const;

    /**
     * @brief Checks if the data of this LAMAArray is available at the passed context.
     *
     * @param[in] context   the context to check for availability of data.
     * @return              if the data of this LAMAArray is available at context
     */
    bool isAvailableAt( ContextPtr context ) const;

    /**
     * @brief Gets the first context where the data of this LAMAArray is available.
     *
     * If possible a context of the passed preferred type is returned.
     *
     * @param[in] preferredType the preferred type for the valid context.
     * @return                  a context there the data of this LAMAArray is available.
     */
    ContextPtr getValidContext( const Context::ContextType preferredType = Context::Host ) const;

    /**
     * Waits for potentially running prefetch.
     */
    void wait() const;

    /**
     * @brief sets the size of this to 0 but does not free any memory
     *
     * This operation should be used before a write-only access that is
     * followed by a resize of the array. It avoids data transfer between
     * two contextes.
     */
    void clear();

    /**
     * @brief sets the size of this to 0 an frees all memory
     */
    void purge();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method.
     */
    virtual Scalar::ScalarType getValueType() const;

    /**
     * @brief reserve a certain amount of data at a specific context
     *
     * @param[in] context where a certain amount of data should be reserved
     * @param[in] capacity amount of data to be allocated
     *
     */
    void reserve( ContextPtr context, const IndexType capacity );

protected:

    using _LAMAArray::mSize;
    using _LAMAArray::constFlag;

    ValueType* get( size_t index );

    const ValueType* get( const size_t index ) const;

    void clear( const size_t index );

    void resize( const size_t index, const IndexType newSize );

    void reserve( const size_t index, const IndexType capacity, const bool copy ) const;

    IndexType capacity( const size_t index ) const;

    /** Complete handling to get read access for a certain context.
     *
     * \returns index of context data array that contains the valid entry.
     */
    int acquireReadAccess( ContextPtr context ) const;

    void releaseReadAccess( const size_t index ) const;

    /** Complete handling to get write access for a certain context.
     *
     * \returns index of context data array that contains the valid entry.
     */
    int acquireWriteAccess( ContextPtr context, bool keepFlag );

    int acquireWriteAccess();

    void releaseWriteAccess( const size_t index );

    void copy( const int toIndex, const int fromIndex ) const;

    /** Copy data from a valid source context to an invalid target context.
     *
     *  This routine assumes that there is already sufficient memory allocated
     *  in the target context.
     */

    void fetch( Context::ContextData& target, const Context::ContextData& source ) const;

    /** Asynchronous copy data from a valid source context to an invalid target context. */

    SyncToken* fetchAsync( Context::ContextData& target, const Context::ContextData& source ) const;

    void getAccess(
        size_t& contextIndex,
        size_t& validIndex,
        ContextPtr context,
        Context::ContextData::AccessKind kind ) const;

    void setHostContext();

    mutable std::vector<Context::ContextData*> mContextData; // available context, pointers are never NULL

    mutable std::auto_ptr<SyncToken> mSyncToken; //!<  outstanding transfers

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

mutable    boost::recursive_mutex access_mutex; // needed to make accesses thread-safe
};

/**
 * @brief LAMAArrayRef is a container that uses external data.
 *
 * @tparam ValueType is the type stored in this container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArrayRef: public LAMAArray<ValueType>
{
public:

    /** Contruct a container for a host array. */

    LAMAArrayRef( ValueType* pointer, IndexType size );

    /** Contruct a container for a const host array.
     *  Due to the const pointer it is guaranteed that the array cannot be modified
     */

    LAMAArrayRef( const ValueType* pointer, IndexType size );

protected:

    using LAMAArray<ValueType>::mSize;

    using LAMAArray<ValueType>::mContextData;
    using LAMAArray<ValueType>::constFlag;
};

/* ---------------------------------------------------------------------------------*/

inline IndexType _LAMAArray::size() const
{
    return mSize;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
template<typename OtherValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const OtherValueType* const values )
                : _LAMAArray( n ), mSyncToken( 0 )
{
    setHostContext();

    if( n <= 0 )
    {
        LAMA_LOG_DEBUG( logger, "Zero-sized array with value constructed: " << *this )
        return;
    }

    Context::ContextData& host = *mContextData[0];
    host.allocate( mSize * sizeof(ValueType) );
    LAMA_LOG_DEBUG( logger, "constructed: " << *this )

    ValueType* host_pointer = static_cast<ValueType*>( host.pointer );

#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( int i = 0; i < mSize; ++i )
    {
        host_pointer[i] = static_cast<ValueType>( values[i] );
    }

    host.valid = true;
    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

}

#endif // LAMA_LAMA_ARRAY_HPP_
