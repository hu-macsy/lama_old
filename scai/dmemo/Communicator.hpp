/**
 * @file Communicator.hpp
 *
 * @license
 * Copyright (c) 2009-2014
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
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/Printable.hpp>

#include <scai/dmemo/CommunicationPlan.hpp>

// internal scai libraris
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/preprocessor.hpp>

// std
#include <memory>
#include <vector>

//#include <cmath>

namespace scai
{

namespace hmemo
{
template<typename ValueType> class HArray;
}

namespace dmemo
{

// Forward declaration of all classes that are used in the interface

class Distribution;

class Halo;

class Communicator;

typedef common::shared_ptr<const Communicator> CommunicatorPtr;

namespace communicator
{

typedef enum
{
    NO,
    MPI,
    GPI,
    MAX_COMMUNICATOR
} CommunicatorKind;

}
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
class COMMON_DLL_IMPORTEXPORT Communicator:

    public  scai::common::Printable,
    public  scai::common::Factory<communicator::CommunicatorKind, CommunicatorPtr>,
    private scai::common::NonCopyable
{

public:

    /** Get a communicator of a certain type from the factory.
     *
     *  @param type is the name of the needed communicator.
     *  @returns pointer to the desired communicator, the default one if not found
     *
     *  More convenient than Factory::create that throws Exception
     */

    static CommunicatorPtr getCommunicator( const communicator::CommunicatorKind& type );

    /** Get a default communicator from the factory.
     *
     *  @returns shared pointer to the default communicator.
     *
     *  The rules for choosing the default communicator are as follows:
     *   - take the one specified by environment/argument SCAI_COMM
     *   - take if available in this order: MPI, GPI, NO
     */

    static CommunicatorPtr getCommunicator();

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

    /**************************************************************************************
     *                                                                                    *
     *  Virtual routines provided by derived class for each supported TypeID              *
     *                                                                                    *
     *************************************************************************************/

    /* @brief Exchange of data between all processors by communication plans.
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
     *
     *  virtual void exchangeByPlan(
     *      TypeId recvData[],
     *      const CommunicationPlan& recvPlan,
     *      const TypeId sendData[],
     *      const CommunicationPlan& sendPlan ) const = 0;
     */

    /* @brief Asynchronous exchange of data between all processors by communication plans.
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
     *
     *  virtual SyncToken* exchangeByPlanAsync(
     *      TypeId recvData[],
     *      const CommunicationPlan& recvPlan,
     *      const TypeId sendData[],
     *      const CommunicationPlan& sendPlan ) const = 0;
     */

    /* @brief Broadcast a typed array from root to all other processors.
     *
     *  @param[in,out] val  in on root, out on all other processors
     *  @param[in]     n    number of elements in vector val
     *  @param[in]     root processor with correct values of val
     *
     *  virtual void bcast( TypeId val[], const IndexType n, const PartitionId root ) const = 0;
     */

    /* @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myvals values that I receive
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allvals values for all processors (size must be n * size() )
     *
     *  virtual void scatter(
     *      TypeId myvals[],
     *      const IndexType n,
     *      const PartitionId root,
     *      const TypeId allvals[] ) const = 0;
     */

    /* @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myvals values that I receive
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allvals values for all processors (size must be sum(sizes) )
     *  @param[in]    sizes   number of total values for all processors
     *
     *  virtual void scatterV(
     *      int myvals[],
     *      const IndexType n,
     *      const PartitionId root,
     *      const int allvals[],
     *      const IndexType sizes[] ) const = 0;
     */

    /* @brief Gather of an array of values from all processors to root.
     *
     *  @param[out]   allvals values that I receive (size must be n * size() )
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myvals values that this processor contributes
     *
     *  virtual void gather(
     *      TypeId allvals[],
     *      const IndexType n,
     *      const PartitionId root,
     *      const TypeId myvals[] ) const = 0;
     */

    /* @brief Gather of an array of double values from all processors to root.
     *
     *  @param[out]   allvals values that I receive (size must be sum(sizes) )
     *  @param[in]    n      number of elements in myvals
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myvals values on this processor
     *  @param[in]    sizes  array with sizes from each processor
     *
     *  The array sizes must only be valid on the root processor. The
     *  output array allvals will only be valid on the root processor.
     *
     *  virtual void gatherV(
     *      TypeId allvals[],
     *      const IndexType n,
     *      const PartitionId root,
     *      const TypeId myvals[],
     *      const IndexType sizes[] ) const = 0;
     */

    /* @brief Sum operations sum up one single value from each partition to a global value.
     *
     *  @param[in] value  value on the calling partition
     *  @returns   global value, available for all partitions.
     *
     *  virtual TypeId sum( const TypeId value ) const = 0;
     *  virtual TypeId min( const TypeId value ) const = 0;
     *  virtual TypeId max( const TypeId value ) const = 0;
     */

    /* @brief Maximal value combined with a location value where maximum was found.
     *
     *  @param[in,out] val        is a value on each processor, only out for root with maximal value
     *  @param[in,out] location   is an additional int value, only out for root
     *  @param[in]     root       rank of processor that has valid results at end of the routine
     *
     *  Only root processor will contain the maximal value and the location loc,
     *
     * virtual void maxloc( TypeId& val, IndexType& location, const PartitionId root ) const = 0;
     */

    /* @brief Swap of an array with another processor.
     *
     * @param[in,out] val is the data array to be swapped
     * @param[in] n is the number of entries in array val
     * @param[in] partner is the rank of partition with which this partition swaps
     *
     * This method can also be used if partner is same as this processor.
     *
     * virtual void swap( TypeId val[], const IndexType n, const PartitionId partner ) const = 0;
     */

    /* @brief This routine shifts data between neighbored processors.
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
     *  virtual IndexType shiftData(
     *      TypeId newVals[],
     *      const IndexType newSize,
     *      const TypeId oldVals[],
     *      const IndexType oldSize,
     *      const int direction ) const = 0;
     */

    /* @brief Asynchronous version of shift.
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
     *
     *  virtual SyncToken* shiftDataAsync(
     *      TypeId newVals[],
     *      const TypeId oldVals[],
     *      const IndexType size,
     *      const int direction ) const = 0;
     */

    /* @brief Exchange of data between all processors by communication plans.
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

    /** All-to-all exchange of an IndexType value between all processors.
     *
     * @param[out] recvValues   will contain one value from each processor
     * @param[in]  sendValues   must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    virtual void all2all( IndexType recvValues[], const IndexType sendValues[] ) const = 0;

    /** @brief Expanded macro by BOOST_PP_REPEAT that defines all virtual routines for one type.
     *
     *  @param z next available repetition dimension
     *  @param I is the actual instantation by BOOST_PP_REPEAT
     *  @param _ remains unused here
     *
     *  Note: this macro is needed as virtual routines must not be templates
     */

#define COMMUNICATOR_METHODS( z, I, _)                                    \
    \
    virtual void exchangeByPlan(                                          \
            ARRAY_TYPE##I* const recvData,                                    \
            const CommunicationPlan& recvPlan,                                \
            const ARRAY_TYPE##I* const sendData,                              \
            const CommunicationPlan& sendPlan ) const = 0;                    \
    \
    virtual tasking::SyncToken* exchangeByPlanAsync(                          \
            ARRAY_TYPE##I* const recvData,                                    \
            const CommunicationPlan& recvPlan,                                \
            const ARRAY_TYPE##I* const sendData,                              \
            const CommunicationPlan& sendPlan ) const = 0;                    \
    \
    virtual void bcast(                                                   \
            ARRAY_TYPE##I val[],                                              \
            const IndexType n,                                                \
            const PartitionId root ) const = 0;                               \
    \
    virtual void all2allv(                                                \
            ARRAY_TYPE##I* recvVal[],                                         \
            IndexType recvCount[],                                            \
            ARRAY_TYPE##I* sendVal[],                                        \
            IndexType sendCount[] ) const = 0;                                \
    \
    virtual void scatter(                                                 \
            ARRAY_TYPE##I myvals[],                                           \
            const IndexType n,                                                \
            const PartitionId root,                                           \
            const ARRAY_TYPE##I allvals[] ) const = 0;                        \
    \
    virtual void scatterV(                                                \
            ARRAY_TYPE##I myvals[],                                           \
            const IndexType n,                                                \
            const PartitionId root,                                           \
            const ARRAY_TYPE##I allvals[],                                    \
            const IndexType sizes[] ) const = 0;                              \
    \
    virtual void gather(                                                  \
            ARRAY_TYPE##I allvals[],                                          \
            const IndexType n,                                                \
            const PartitionId root,                                           \
            const ARRAY_TYPE##I myvals[] ) const = 0;                         \
    \
    virtual void gatherV(                                                 \
            ARRAY_TYPE##I allvals[],                                          \
            const IndexType n,                                                \
            const PartitionId root,                                           \
            const ARRAY_TYPE##I myvals[],                                     \
            const IndexType sizes[] ) const = 0;                              \
    \
    virtual void swap(                                                    \
            ARRAY_TYPE##I val[],                                              \
            const IndexType n,                                                \
            const PartitionId partner ) const = 0;                            \
    \
    virtual void maxloc(                                                  \
            ARRAY_TYPE##I& val,                                               \
            IndexType& location,                                              \
            const PartitionId root ) const = 0;                               \
    \
    virtual ARRAY_TYPE##I min(                                            \
            const ARRAY_TYPE##I value ) const = 0;                            \
    \
    virtual ARRAY_TYPE##I sum(                                            \
            const ARRAY_TYPE##I value ) const = 0;                            \
    \
    virtual ARRAY_TYPE##I max(                                            \
            const ARRAY_TYPE##I value ) const = 0;                            \
    \
    virtual IndexType shiftData(                                          \
            ARRAY_TYPE##I newVals[],                                          \
            const IndexType newSize,                                          \
            const ARRAY_TYPE##I oldVals[],                                    \
            const IndexType oldSize,                                          \
            const int direction ) const = 0;                                  \
    \
    virtual tasking::SyncToken* shiftDataAsync(                               \
            ARRAY_TYPE##I newVals[],                                          \
            const ARRAY_TYPE##I oldVals[],                                    \
            const IndexType size,                                             \
            const int direction ) const = 0;                                  \
     

    // define communicator methods for all supported data types

    BOOST_PP_REPEAT( ARRAY_TYPE_CNT, COMMUNICATOR_METHODS, _ )

#undef COMMUNICATOR_METHODS

    /**************************************************************************************
     *                                                                                    *
     *  Other Virtual routines provided by derived class                                  *
     *                                                                                    *
     *************************************************************************************/

    // Sum of size_t values needed for large quantities like memory usage
    virtual size_t sum( const size_t value ) const = 0;

    // Broadcast of characters needed for strings

    virtual void bcast( char val[], const IndexType n, const PartitionId root ) const = 0;

    /** @brief Barrier synchronization between all processors. */

    virtual void synchronize() const = 0;

    /**************************************************************************************
     *                                                                                    *
     *  Communication routines provided by this base class using other routines           *
     *                                                                                    *
     *************************************************************************************/

    // Boolean global reduction operations are implemented by using the sum reduction.
    virtual bool all( const bool flag ) const;
    virtual bool any( const bool flag ) const;

    /** Broadcast of a string.
     *
     *  This class provides one implementation, but derived classes might override it.
     */
    virtual void bcast( std::string&, const PartitionId root ) const;


    template<typename ValueType>
    void all2allv( ValueType* recvBuffer[], IndexType recvCount[],  ValueType* sendBuffer[], IndexType sendCount[] ) const;

    /** @brief allgather is combination of gather and broadcast
     *
     *  @param[out]  allvals result array with values from all processors ( allocated for at least n * getSize() )
     *  @param[in]   n      number of elements in myvals
     *  @param[in]   myvals array with n values that this processor contributes
     *
     *  This routine will be implemented by other available routines.
     */
    template<typename ValueType>
    void allgather( ValueType allvals[], const IndexType n, const ValueType myvals[] ) const
    {
        PartitionId root = 0;
        gather( allvals, n, root, myvals );
        bcast( allvals, n * getSize(), root );

        // @ToDo: implement this by a circular shift
    }

    /** Ownership computation for indexes where only each partition
     * individually can determine whether an index is local or not.
     *
     * This method should not be called for distributions where the
     * owner can be computed by a closed formula.
     */
    void computeOwners(
        PartitionId owners[],
        const Distribution& distribution,
        const IndexType requiredIndexes[],
        const IndexType n ) const;

    void computeOwners(
        const std::vector<IndexType>& requiredIndexes,
        const Distribution& distribution,
        std::vector<PartitionId>& owners ) const
    {
        IndexType n = requiredIndexes.size();
        owners.resize( n );
        computeOwners( &owners[0], distribution, &requiredIndexes[0], n );
    }

    /**************************************************************************************
     *                                                                                    *
     *  Communication routines for HArrays for convenience                             *
     *                                                                                    *
     *************************************************************************************/

    /** exchangeByPlan for HArrays instead of usual array */

    template<typename ValueType>
    void exchangeByPlan(
        hmemo::HArray<ValueType>& recvArray,
        const CommunicationPlan& recvPlan,
        const hmemo::HArray<ValueType>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** Asynchronous exchange of HArrays. */

    template<typename ValueType>
    tasking::SyncToken* exchangeByPlanAsync(
        hmemo::HArray<ValueType>& recvArray,
        const CommunicationPlan& recvPlan,
        const hmemo::HArray<ValueType>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** @brief Update of halo array via Halo object.
     *
     *  @tparam     ValueType             arithmetic type of involved arrays
     *  @param[out] haloValues    will contain the non-local values from other processors
     *  @param[in]  localValues   is the local part of the array on each processor
     *  @param[in]  halo is the   Halo object containing all information about exchange
     *
     *  This method is not virtual but will use the pure virtual methods of base classes.
     */
    template<typename ValueType>
    void updateHalo(
        hmemo::HArray<ValueType>& haloValues,
        const hmemo::HArray<ValueType>& localValues,
        const Halo& halo ) const;

    /** @brief Asynchronous update of halo array via Halo object. */

    template<typename ValueType>
    tasking::SyncToken* updateHaloAsync(
        hmemo::HArray<ValueType>& haloValues,
        const hmemo::HArray<ValueType>& localValues,
        const Halo& halo ) const;

    /** @brief Shift on LAMA arrays.
     *
     *  @tparam     ValueType           type of data to be shifted
     *  @param[out] recv        array to receive  this partition
     *  @param[in]  send        array to send from this partition
     *  @param[in]  direction   number of positions to shift, e.g. 1 or -1
     *
     *  Note: The recv array must have a capacity that is sufficent to
     *        receive all the data.
     */
    template<typename ValueType>
    void shiftArray( hmemo::HArray<ValueType>& recv, const hmemo::HArray<ValueType>& send, const int direction ) const;

    /** @brief Asychronous shift on LAMA arrays.
     *
     *  @tparam     ValueType           TODO[doxy] Complete Description.
     *  @param[out] recvArray   array to receive for this partition
     *  @param[in]  sendArray   array to send from this partition
     *  @param[in]  direction   number of positions to shift, e.g. 1 or -1
     *
     *  Note: All partitions must have the same size for send/recv array
     */
    template<typename ValueType>
    tasking::SyncToken* shiftAsync(
        hmemo::HArray<ValueType>& recvArray,
        const hmemo::HArray<ValueType>& sendArray,
        const int direction ) const;

    /** Override routine of base class Printable. */

    virtual void writeAt( std::ostream& stream ) const;

    /** @brief Getter for the type of a communicator. */

    const communicator::CommunicatorKind& getType() const
    {
        return mCommunicatorType;
    }

    /** This routine provides a context that might be most efficient for
     *  the communication of data.
     *
     * @param[in] array is a LAMA array that is used for communication
     * @returns a context pointer where data should be available for communication
     *
     * Note: this routine is mainly used for sending of values; in case of CUDA aware
     *       communication it might return a CUDA context if the array has valid data there.
     */
    virtual hmemo::ContextPtr getCommunicationContext( const hmemo::_HArray& array ) const = 0;

protected:

    // Default constructor can only be called by base classes.

    Communicator( const communicator::CommunicatorKind& type );

    communicator::CommunicatorKind mCommunicatorType;

    int mNodeRank; // rank of this processor on its node
    int mNodeSize; // number of processors on same node

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Read in the environment variable LAMA_NP for user processor array.
     *
     * @param[out] userProcArray specifies the user processor array.
     *
     */
    static    void getUserProcArray( PartitionId userProcArray[3] );

    /** Shift implementation for direction == 0, just copies values. */

    template<typename ValueType>
    IndexType shift0( ValueType newVals[], const IndexType newSize,
                      const ValueType oldVals[], const IndexType oldSize ) const;

};

/* -------------------------------------------------------------------------- */

PartitionId Communicator::getNeighbor( int pos ) const
{
    PartitionId size = getSize();
    PartitionId rank = getRank();

    PartitionId apos = common::Math::abs( pos );

    SCAI_ASSERT( apos <= size, "neighbor pos " << pos << " out of range (" << size << ")" )

    return ( size + rank + pos ) % size;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::exchangeByPlan(
    hmemo::HArray<ValueType>& recvArray,
    const CommunicationPlan& recvPlan,
    const hmemo::HArray<ValueType>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    SCAI_ASSERT_EQ_ERROR( sendArray.size(), sendPlan.totalQuantity(), "size mismatch" )

    IndexType recvSize = recvPlan.totalQuantity();

    // find a context where data of sendArray can be communicated
    // if possible try to find a context where valid data is available
    // CUDAaware MPI: might give GPU or Host context here

    hmemo::ContextPtr comCtx = getCommunicationContext( sendArray );

    SCAI_LOG_DEBUG( logger, *this << ": exchangeByPlan, comCtx = " << *comCtx )

    hmemo::ReadAccess<ValueType> sendData( sendArray, comCtx );

    // Data will be received at the same context where send data is

    hmemo::WriteOnlyAccess<ValueType> recvData( recvArray, comCtx, recvSize );

    exchangeByPlan( recvData.get(), recvPlan, sendData.get(), sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* Communicator::exchangeByPlanAsync(
    hmemo::HArray<ValueType>& recvArray,
    const CommunicationPlan& recvPlan,
    const hmemo::HArray<ValueType>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    SCAI_ASSERT_EQ_ERROR( sendArray.size(), sendPlan.totalQuantity(), "size mismatch" )

    IndexType recvSize = recvPlan.totalQuantity();

    // allocate accesses, SyncToken will take ownership

    hmemo::ContextPtr comCtx = getCommunicationContext( sendArray );

    SCAI_LOG_DEBUG( logger, *this << ": exchangeByPlanAsync, comCtx = " << *comCtx )

    hmemo::ReadAccess<ValueType> sendData( sendArray, comCtx );
    hmemo::WriteOnlyAccess<ValueType> recvData( recvArray, comCtx, recvSize );

    tasking::SyncToken* token( exchangeByPlanAsync( recvData.get(), recvPlan, sendData.get(), sendPlan ) );

    // Add the read and write access to the sync token to get it freed after successful wait
    // conversion common::shared_ptr<hmemo::HostWriteAccess<ValueType> > -> common::shared_ptr<BaseAccess> supported

    token->pushRoutine( recvData.releaseDelayed() );
    token->pushRoutine( sendData.releaseDelayed() );

    // return ownership of new created object

    return token;
}

} /* end namespace dmemo */

} /* end namespace scai */
