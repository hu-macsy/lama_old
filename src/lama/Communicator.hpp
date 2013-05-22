/**
 * @file Communicator.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Base and interface class for communicators used in LAMA
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 * $Id$
 */
#ifndef LAMA_COMMUNICATOR_HPP_
#define LAMA_COMMUNICATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>
#include <lama/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/CommunicationPlan.hpp>
#include <lama/SyncToken.hpp>

#include <lama/exception/LAMAAssert.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

#include <memory>
#include <vector>
#include <cmath>

namespace lama
{

// Forward declaration of all classes that are used in the interface

template<typename T> class LAMAArray;

class Distribution;

class Halo;

/**
 * @brief Base and interface class for communicators used in LAMA
 *
 * This base class defines pure methods for communication and data exchange
 * between the different partitions of a distributed application. These pure
 * methods must be provided by base classes that implement a full communicator.
 *
 * Furthermore, this class provides also some higher functionality communication methods
 * especially for LAMA arrays that will use the pure virtual methods.
 *
 * Communicator objects will be managed by the Communicator factory. All
 * operations on a communicator are const methods. Communicators are referenced
 * via shared pointers so that resources are freed only if no more reference
 * to the communicator exists.
 *
 * Default copy constructor and assignment operator are disabled.
 *
 */
class LAMA_DLL_IMPORTEXPORT Communicator: public Printable, private NonCopyable
{

public:

    /** Enumeration type for supported thread safety levels. */

    enum ThreadSafetyLevel
    {
        Funneled = 1, //!< The Program can be multithreaded but only the master thread should communicate
        Serialized = 2, //!< The Program can be multithreaded but the communicator usage is serialized
        Multiple = 3 //!< The Program can be multithreaded and the communicator can be used by all threads concurrently
    };

    virtual ~Communicator();

    /** @brief TODO[doxy] Complete Description.
     *
     *  @param[in] sizeX      TODO[doxy] Complete Description.
     *  @param[in] sizeY      TODO[doxy] Complete Description.
     *  @param[in] procgrid   TODO[doxy] Complete Description.
     */
    void factorize2( const double sizeX, const double sizeY, PartitionId procgrid[2] ) const;

    /** @brief TODO[doxy] Complete Description.
     *
     *  @param[in] sizeX      TODO[doxy] Complete Description.
     *  @param[in] sizeY      TODO[doxy] Complete Description.
     *  @param[in] sizeZ      TODO[doxy] Complete Description.
     *  @param[in] procgrid   TODO[doxy] Complete Description.
     */
    void factorize3( const double sizeX, const double sizeY, const double sizeZ, PartitionId procgrid[3] ) const;

    /** @brief TODO[doxy] Complete Description.
     *
     *  @param[out] pos       TODO[doxy] Complete Description.
     *  @param[in]  procgrid  TODO[doxy] Complete Description.
     */
    void getGrid2Rank( PartitionId pos[2], const PartitionId procgrid[2] ) const;

    /** @brief TODO[doxy] Complete Description.
     *
     *  @param[out] pos       TODO[doxy] Complete Description.
     *  @param[in]  procgrid  TODO[doxy] Complete Description.
     */
    void getGrid3Rank( PartitionId pos[3], const PartitionId procgrid[3] ) const;

    /** @brief Equality operator for two communicators.
     *
     *  @param[in] other   communicator for comparison
     *  @return            true if this communicator is equal to other communicator.
     *
     *  Note: Logic of this operator is implemented via the virtual function
     *  isEqual.
     */
    bool operator==( const Communicator& other ) const;

    /** @brief Diversity operator for two communicators.
     *
     *  @param[in] other   communicator for comparison
     *  @returns           true if this communicator is unequal to other communicator.
     *
     *  Note: Logic of this operator is implemented via the virtual function
     *  isEqual.
     */
    bool operator!=( const Communicator& other ) const;

    /** @brief Virtual method used for the equality operator. */

    virtual bool isEqual( const Communicator& other ) const = 0;

    /** @brief Getter function for the thread safety level.
     *
     *  @return   the thread safety level.
     */

    virtual ThreadSafetyLevel getThreadSafetyLevel() const = 0;

    /** @brief Getter of the number of partitions.
     *
     * @return number of partitions.
     */

    virtual PartitionId getSize() const = 0;

    /** @brief Getter of rank of this partition.
     *
     * @return rank of this partition, 0 <= rank < getSize()
     */
    virtual PartitionId getRank() const = 0;

    /** @brief Getter of the number of partitions on same node.
     *
     * @return number of partitions on same node as this partition.
     */

    virtual PartitionId getNodeSize() const = 0;

    /** @brief Getter of node rank 
     *
     * @return rank of this partition on its node, 0 <= rank < getNodeSize()
     */
    virtual PartitionId getNodeRank() const = 0;

    /**
     * Help routine to get the rank of a neighbored position.
     *
     * @param[in] pos   is the distance to the neighbor (also negative)
     * @return          rank of the neighbored processor
     *
     * This method assumes a circular ring formed by all processors.
     */
    inline PartitionId getNeighbor( int pos ) const;

    /** All-to-all exchange of an integer value between all processors.
     *
     * @param[out] recvValues   will contain one value from each processor
     * @param[in]  sendValues   must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( int* recvValues, const int* sendValues ) const = 0;

    /** @brief Exchange of data between all processors by communication plans.
     *
     *  @param[out] recvData   buffer for data received from other processors
     *  @param[in]  recvPlan   number of elements and offsets for receiving
     *  @param[in]  sendData   buffer for data to send to other processors
     *  @param[in]  sendPlan   contains number of elements and offsets for sending data
     *
     *  All send and receive data between each pair of processors must be a contiguous
     *  part of the sendData or recvData.
     *
     *  The size of recvData must be recvPlan.totalQuantity().
     *  The size of sendData must be sendPlan.totalQuantity().
     */
    virtual void exchangeByPlan(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Exchange of float values. */

    virtual void exchangeByPlan(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Exchange of double values. */

    virtual void exchangeByPlan(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Exchange of LAMAArrays. */

    template<typename T>
    void exchangeByPlan(
        LAMAArray<T>& recvArray,
        const CommunicationPlan& recvPlan,
        const LAMAArray<T>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** @brief Asynchronous exchange of data between all processors by communication plans.
     *
     *  @param[out] recvData   buffer for data received from other processors
     *  @param[in]  recvPlan   number of elements and offsets for receiving
     *  @param[in]  sendData   buffer for data to send to other processors
     *  @param[in]  sendPlan   contains number of elements and offsets for sending data
     *
     *  @return TODO[doxy] Complete Description.
     *
     *  All send and receive data between each pair of processors must be a contiguous
     *  part of the sendData or recvData.
     *
     *  The size of recvData must be recvPlan.totalQuantity().
     *  The size of sendData must be sendPlan.totalQuantity().
     */
    virtual SyncToken* exchangeByPlanAsync(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Asynchronous exchange of float values. */

    virtual SyncToken* exchangeByPlanAsync(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Asynchronous exchange of double values. */

    virtual SyncToken* exchangeByPlanAsync(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const = 0;

    /** Asynchronous exchange of LAMAArrays. */

    template<typename T>
    SyncToken* exchangeByPlanAsync(
        LAMAArray<T>& recvArray,
        const CommunicationPlan& recvPlan,
        const LAMAArray<T>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** @brief Update of halo array via Halo object.
     *
     *  @tparam     T             TODO[doxy] Complete Description.
     *  @param[out] haloValues    will contain the non-local values from other processors
     *  @param[in]  localValues   is the local part of the array on each processor
     *  @param[in]  halo is the   Halo object containing all information about exchange
     *
     *  This method is not virtual but will use the pure virtual methods of base classes.
     */
    template<typename T>
    void updateHalo( LAMAArray<T>& haloValues, const LAMAArray<T>& localValues, const Halo& halo ) const;

    /** @brief Asynchronous update of halo array via Halo object. */

    template<typename T>
    SyncToken* updateHaloAsync(
        LAMAArray<T>& haloValues,
        const LAMAArray<T>& localValues,
        const Halo& halo ) const;

    /** @brief Shift on LAMA arrays.
     *
     *  @tparam     T           type of data to be shifted       
     *  @param[out] recv        array to receive  this partition
     *  @param[in]  send        array to send from this partition
     *  @param[in]  direction   number of positions to shift, e.g. 1 or -1
     *
     *  Note: The recv array must have a capacity that is sufficent to
     *        receive all the data.
     */
    template<typename T>
    void shiftArray( LAMAArray<T>& recv, const LAMAArray<T>& send, const int direction ) const;

    /** @brief Asychronous shift on LAMA arrays.
     *
     *  @tparam     T           TODO[doxy] Complete Description.
     *  @param[out] recvArray   array to receive for this partition
     *  @param[in]  sendArray   array to send from this partition
     *  @param[in]  direction   number of positions to shift, e.g. 1 or -1
     *
     *  Note: All partitions must have the same size for send/recv array
     */
    template<typename T>
    SyncToken* shiftAsync(
        LAMAArray<T>& recvArray,
        const LAMAArray<T>& sendArray,
        const int direction ) const;

    /** Ownership computation for indexes where only each partition
     * individually can determine whether an index is local or not.
     *
     * This method should not be called for distributions where the
     * owner can be computed by a closed formula.
     */
    void computeOwners(
        const std::vector<IndexType>& requiredIndexes,
        const Distribution& distribution,
        std::vector<PartitionId>& owners ) const;

    /** @brief Broadcast a typed array from root to all other processors.
     *
     *  @param[in,out] val  in on root, out on all other processors
     *  @param[in]     n    number of elements in vector val
     *  @param[in]     root processor with correct values of val
     */
    virtual void bcast( double val[], const IndexType n, const PartitionId root ) const = 0;
    virtual void bcast( float val[], const IndexType n, const PartitionId root ) const = 0;
    virtual void bcast( int val[], const IndexType n, const PartitionId root ) const = 0;

    /** @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myvals values that I receive
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allvals values for all processors (size must be n * size() )
     */
    virtual void scatter(
        double myvals[],
        const IndexType n,
        const PartitionId root,
        const double allvals[] ) const = 0;
    virtual void scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const = 0;
    virtual void scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const = 0;

    /** @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myvals values that I receive
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allvals values for all processors (size must be sum(sizes) )
     *  @param[in]    sizes   number of total values for all processors
     */
    virtual void scatter(
        double myvals[],
        const IndexType n,
        const PartitionId root,
        const double allvals[],
        const IndexType sizes[] ) const = 0;
    virtual void scatter(
        float myvals[],
        const IndexType n,
        const PartitionId root,
        const float allvals[],
        const IndexType sizes[] ) const = 0;
    virtual void scatter(
        int myvals[],
        const IndexType n,
        const PartitionId root,
        const int allvals[],
        const IndexType sizes[] ) const = 0;

    /** @brief Gather of an array of values from all processors to root.
     *
     *  @param[out]   allvals values that I receive (size must be n * size() )
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myvals values for all processors
     */
    virtual void gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const = 0;
    virtual void gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const = 0;
    virtual void gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const = 0;

    /** @brief Gather of an array of values from all processors to root.
     *
     *  @param[out]   allvals values that I receive (size must be sum(sizes) )
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myvals values for all processors
     *  @param[in]    sizes   number of total values for all processors
     */
    virtual void gather(
        double allvals[],
        const IndexType n,
        const PartitionId root,
        const double myvals[],
        const IndexType sizes[] ) const = 0;
    virtual void gather(
        float allvals[],
        const IndexType n,
        const PartitionId root,
        const float myvals[],
        const IndexType sizes[] ) const = 0;
    virtual void gather(
        int allvals[],
        const IndexType n,
        const PartitionId root,
        const int myvals[],
        const IndexType sizes[] ) const = 0;

    /** @brief This routine shifts data between neighbored processors.
     *
     *  @param[out] newVals  array with data this partition get from neighbored partition
     *  @param[in]  newSize  allocated size of array newVals
     *  @param[in]  oldVals  array with data this partition sends to neighbored partition
     *  @param[in]  oldSize  number of elements of array oldVals that will be sent
     *  @param[in]  direction specifies the neighbored partitions to send and receive
     *
     *  @returns number of values received in newVals (newSize must have been larger before)
     *
     *  Each partition sends to rank() + direction and receives from rank() - direction.
     *
     */
    template <typename T>
    IndexType shiftData( T newVals[], const IndexType newSize,
                         const T oldVals[], const IndexType oldSize,
                         const int direction ) const
    {
        return shiftImpl( newVals, newSize, oldVals, oldSize, direction );
    }
 
    /** @brief Asynchronous version of shift.
     *
     *  @param[out] newVals  array with data this partition get from neighbored partition
     *  @param[in]  oldVals  array with data this partition sends to neighbored partition
     *  @param[in]  size  number of elements of array oldVals to be sent andd newVals to receive
     *  @param[in]  direction specifies the neighbored partitions to send and receive
     *
     *  As there is no information about the received size this routine can only be called for
     *  arrays that have the same size on all partitions.
     *
     *  A default implementation is provided that returns a NoSyncToken. Derived classes
     *  should override this method if there is a benefit of using asynchronous transfers.
     */

    template <typename T>
    SyncToken* shiftDataAsync(
        T newVals[],
        const T oldVals[],
        const IndexType size,
        const int direction ) const
    {
        // call virtual implementation routine 

        return shiftAsyncImpl( newVals, oldVals, size, direction );
    }

    /** @brief Sum operations sum up one single value from each partition to a global value.
     *
     *  @param[in] value  value on the calling partition
     *  @returns   global value, available for all partitions.
     */
    virtual float sum( const float value ) const = 0;
    virtual double sum( const double value ) const = 0;
    virtual int sum( const int value ) const = 0;
    virtual size_t sum( const size_t value ) const = 0;

    virtual float min( const float value ) const = 0;
    virtual float max( const float value ) const = 0;

    virtual double min( const double value ) const = 0;
    virtual double max( const double value ) const = 0;

    virtual int min( const int value ) const = 0;
    virtual int max( const int value ) const = 0;

    // Boolean reduction operations can be implemented by using the sum reduction.

    virtual bool all( const bool flag ) const;
    virtual bool any( const bool flag ) const;

    /** @brief Maximal value combined with a location value where maximum was found.
     *
     *  @param[in,out] val        is a value on each processor, only out for root with maximal value
     *  @param[in,out] location   is an additional int value, only out for root
     *  @param[in]     root       TODO[doxy] Complete Description.
     *
     *  Only root processor will contain the maximal value and the location loc.
     */
    virtual void maxloc( double& val, int& location, const PartitionId root ) const = 0;
    virtual void maxloc( float& val, int& location, const PartitionId root ) const = 0;
    virtual void maxloc( int& val, int& location, const PartitionId root ) const = 0;

    /** @brief Swap of an array with another processor.
     *
     * @param[in,out] val is the data array to be swapped
     * @param[in] n is the number of entries in array val
     * @param[in] partner is the rank of partition with which this partition swaps
     *
     * This method can also be used if partner is same as this processor.
     */
    virtual void swap( double val[], const IndexType n, const PartitionId partner ) const = 0;
    virtual void swap( float val[], const IndexType n, const PartitionId partner ) const = 0;
    virtual void swap( int val[], const IndexType n, const PartitionId partner ) const = 0;

    /** Gather single float value from each processor into a vector. */

    virtual void gather( std::vector<float>& values, float value ) const = 0;

    /** @brief Barrier synchronization between all processors. */

    virtual void synchronize() const = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /** @brief Getter for the type of a communicator. */

    const std::string& getType() const
    {
        return mCommunicatorType;
    }

protected:

    // Default constructor can only be called by base classes.

    Communicator( const std::string& type );

    std::string mCommunicatorType;

    int mNodeRank; // rank of this processor on its node
    int mNodeSize; // number of processors on same node

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** Read in the environment variable LAMA_NP for user processor array.
     *
     * @param[out] userProcArray specifies the user processor array.
     *
     */
    static void getUserProcArray( PartitionId userProcArray[3] );

    /** Shift implementation for direction == 0, just copies values. */

    template<typename T>
    IndexType shift0( T newVals[], const IndexType newSize,
                      const T oldVals[], const IndexType oldSize ) const;

    /** Default asynchronous shift uses synchronous shift. */

    template<typename T>
    SyncToken* defaultShiftAsync( T newVals[], const T oldVals[],
            const IndexType size, const int direction ) const;

    /** getter for Context needed for Communication
     *
     * @returns the ContextPtr
     */
    virtual ContextPtr getCommunicationContext() const = 0;

    /** Implemenation of shiftData for double, must be provided by derived classes. */

    virtual IndexType shiftImpl(
        double newVals[],
        const IndexType newSize,
        const double oldVals[],
        const IndexType oldSize,
        const int direction ) const = 0;

    /** Implemenation of shiftData for float, must be provided by derived classes. */

    virtual IndexType shiftImpl(
        float newVals[],
        const IndexType newSize,
        const float oldVals[],
        const IndexType oldSize,
        const int direction ) const = 0;

    /** Implemenation of shiftData for int, must be provided by derived classes. */

    virtual IndexType shiftImpl(
        int newVals[],
        const IndexType newSize,
        const int oldVals[],
        const IndexType oldSize,
        const int direction ) const = 0;

    /** Implemenation of shiftDataAsync for int, can be overriden by derived classes. */

    virtual SyncToken* shiftAsyncImpl(
        double newVals[],
        const double oldVals[],
        const IndexType size,
        const int direction ) const;

    /** Implemenation of shiftDataAsync for float, can be overriden by derived classes. */

    virtual SyncToken* shiftAsyncImpl(
        float newVals[],
        const float oldVals[],
        const IndexType size,
        const int direction ) const;

    /** Implemenation of shiftDataAsync for double, can be overriden by derived classes. */

    virtual SyncToken* shiftAsyncImpl(
        int newVals[],
        const int oldVals[],
        const IndexType size,
        const int direction ) const;
};

typedef boost::shared_ptr<const Communicator> CommunicatorPtr;

/* -------------------------------------------------------------------------- */

PartitionId Communicator::getNeighbor( int pos ) const
{
    PartitionId size = getSize();
    PartitionId rank = getRank();

    LAMA_ASSERT( std::abs(pos) <= size, "neighbor pos "<<pos<<" out of range ("<<size<<")" )

    return ( size + rank + pos ) % size;
}

/* -------------------------------------------------------------------------- */

template<typename T>
void Communicator::exchangeByPlan(
    LAMAArray<T>& recvArray,
    const CommunicationPlan& recvPlan,
    const LAMAArray<T>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_ERROR(
        sendArray.size() == sendPlan.totalQuantity(),
        "Send array has size " << sendArray.size() << ", but send plan requires " << sendPlan.totalQuantity() << " entries" )

    IndexType recvSize = recvPlan.totalQuantity();

    ContextPtr comCtx = getCommunicationContext();

    WriteAccess<T> recvData( recvArray, comCtx );
    ReadAccess<T> sendData( sendArray, comCtx );

    recvData.clear();
    recvData.resize( recvSize );

    exchangeByPlan( recvData.get(), recvPlan, sendData.get(), sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename T>
SyncToken* Communicator::exchangeByPlanAsync(
    LAMAArray<T>& recvArray,
    const CommunicationPlan& recvPlan,
    const LAMAArray<T>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_ERROR(
        sendArray.size() == sendPlan.totalQuantity(),
        "Send array has size " << sendArray.size() << ", but send plan requires " << sendPlan.totalQuantity() << " entries" );

    IndexType recvSize = recvPlan.totalQuantity();

    // allocate accesses, SyncToken will take ownership

    boost::shared_ptr<HostWriteAccess<T> > recvData( new HostWriteAccess<T>( recvArray ) );
    boost::shared_ptr<HostReadAccess<T> > sendData( new HostReadAccess<T>( sendArray ) );

    recvData->clear();
    recvData->resize( recvSize );

    SyncToken* token( exchangeByPlanAsync( recvData->get(), recvPlan, sendData->get(), sendPlan ) );

    // Add the read and write access to the sync token to get it freed after successful wait
    // conversion boost::shared_ptr<HostWriteAccess<T> > -> boost::shared_ptr<BaseAccess> supported

    token->pushAccess( recvData ); 
    token->pushAccess( sendData );

    // return ownership of new created object 

    return token;
}

}

#endif // LAMA_COMMUNICATOR_HPP_
