/**
 * @file Communicator.hpp
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
 * @brief Base and interface class for communicators used in LAMA
 * @author Jiri Kraus, Thomas Brandes
 * @date 23.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/dmemo/CommunicationPlan.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/Printable.hpp>

// internal scai libraris
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/count.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/loop.hpp>

// std
#include <memory>
#include <vector>
#include <string>
#include <stack>

//#include <cmath>

namespace scai
{

namespace hmemo
{
template<typename ValueType> class HArray;
}

/** @brief The namespace dmemo contains all stuff belonging to project dmemo (distributed memory). */

namespace dmemo
{

// Forward declaration of all classes that are used in the interface

class Distribution;

class Communicator;

class CommunicationPlan;

typedef std::shared_ptr<const Communicator> CommunicatorPtr;

/** @brief Enum call CommunicatorType and its values */

enum class CommunicatorType
{
    NO,                  //!< No communicator
    MPI,                 //!< MPI communicator
    MAX_COMMUNICATOR     //!< dummy value for number of communicators
};

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const CommunicatorType& kind );

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
    public  scai::common::Factory<CommunicatorType, CommunicatorPtr>,
    private scai::common::NonCopyable,
    public  std::enable_shared_from_this<const Communicator>
{

public:

    /** Get a communicator of a certain type from the factory.
     *
     *  @param type specifies the required communicator.
     *  @returns shared pointer to the queried communicator
     *
     *  @throws common::Exception for unsupported type
     */

    static CommunicatorPtr getCommunicatorPtr( const CommunicatorType& type );

    /** Get communicator from the factory.
     *
     *  @returns shared pointer to the default communicator.
     *
     *  The rules for choosing the default communicator are as follows:
     *   - take the one specified by environment/argument SCAI_COMMUNICATOR
     *   - if not specified take the default one
     */

    static CommunicatorPtr getCommunicatorPtr();

    /** Get the current/actual communicator. */

    static const Communicator& getCurrent();

    /** Get a default communicator from the factory.
     *
     *  @returns shared pointer to the default communicator.
     *
     *  The rules for choosing the default communicator are as follows:
     *   - If MPI is available, use the MPI communicator.
     *   - Otherwise, use no communicator.
     */

    static CommunicatorPtr getDefaultCommunicatorPtr();

    /** Enumeration type for supported thread safety levels. */

    enum ThreadSafetyLevel
    {
        Funneled = 1,   //!< The Program can be multithreaded but only the master thread should communicate
        Serialized = 2, //!< The Program can be multithreaded but the communicator usage is serialized
        Multiple = 3    //!< The Program can be multithreaded and the communicator can be used by all threads concurrently
    };

    virtual ~Communicator();

    /** @brief Find a facotorization of size with two factors
     *
     *  @param[in] size       total number of processors
     *  @param[in] sizeX      weight for the first dimension
     *  @param[in] sizeY      weight for the second dimension
     *  @param[out] procgrid   array of size 2 with the two factors
     */
    static void factorize2( PartitionId procgrid[2], const PartitionId size, const double sizeX, const double sizeY );

    /** @brief Find a factorization of size with three factors
     *
     *  @param[out] procgrid   array of size 3 with the three factors
     *  @param[in] size       total number of processors
     *  @param[in] sizeX      weight for the first dimension
     *  @param[in] sizeY      weight for the second dimension
     *  @param[in] sizeZ      weight for the third dimension
     */
    static void factorize3( PartitionId procgrid[3], const PartitionId size, const double sizeX, const double sizeY, const double sizeZ );

    /** @brief Get the position of this processor in a two-dimensional processor grid
     *
     *  @param[out] pos       my position
     *  @param[in]  procgrid  sizes of the two dimensions
     */
    void getGrid2Rank( PartitionId pos[2], const PartitionId procgrid[2] ) const;

    /** @brief Get the position of this processor in a three-dimensional processor grid
     *
     *  @param[out] pos       my position
     *  @param[in]  procgrid  sizes of the three dimensions
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

    inline PartitionId getSize() const;

    /** @brief Getter of rank of this partition.
     *
     * @return rank of this partition, 0 <= rank < getSize()
     */
    inline PartitionId getRank() const;

    /** @brief Getter of the number of partitions on same node.
     *
     * @return number of partitions on same node as this partition.
     */
    inline PartitionId getNodeSize() const;

    /** @brief Getter of node rank
     *
     * @return rank of this partition on its node, 0 <= rank < getNodeSize()
     */
    inline PartitionId getNodeRank() const;

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

    /** Create new communicators based on colors and keys 
     *
     *  @param[in] color processors with the same color will be in the same communicator
     *  @param[in] key   used to rank processor in new communicator
     *  @return pointer to a new created communicator
     */
    virtual Communicator* splitIt( PartitionId color, PartitionId key ) const = 0;

    CommunicatorPtr split( PartitionId color ) const
    {
        return CommunicatorPtr( splitIt( color, getRank() ) );
    }

    CommunicatorPtr split( PartitionId color, PartitionId key ) const
    {
        return CommunicatorPtr( splitIt( color, key ) );
    }

    /** Build a new communication plan that is the inverse of plan
     *
     *  processor[p].plan->entry[q].quantity  = processor[q].entry[p].quantity
     *
     *  Note: this is just an all to all communication
     */
    CommunicationPlan transpose( const CommunicationPlan& plan ) const; 

    /* @brief Exchange of data between all processors by communication plans.
     *
     *  @param[out] recvVals   buffer for data received from other processors
     *  @param[in]  recvPlan   number of elements and offsets for receiving
     *  @param[in]  sendVals   buffer for data to send to other processors
     *  @param[in]  sendPlan   contains number of elements and offsets for sending data
     *
     *  All send and receive data between each pair of processors must be a contiguous
     *  part of the sendData or recvData.
     *
     *  The size of recvData must be recvPlan.totalQuantity().
     *  The size of sendData must be sendPlan.totalQuantity().
     */
    template<typename ValueType>
    void exchangeByPlan(
        ValueType recvVals[],
        const CommunicationPlan& recvPlan,
        const ValueType sendVals[],
        const CommunicationPlan& sendPlan ) const;

    /** Pure virtual routine, no template arguments as ValueType is coded via ScalarType */

    virtual void exchangeByPlanImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::ScalarType stype ) const = 0;

    /* @brief Asynchronous exchange of data between all processors by communication plans.
     *
     *  @param[out] recvVals   buffer for data received from other processors
     *  @param[in]  recvPlan   number of elements and offsets for receiving
     *  @param[in]  sendVals   buffer for data to send to other processors
     *  @param[in]  sendPlan   contains number of elements and offsets for sending data
     *
     *  @return SyncToken that can be used to wait for completion
     *
     *  All send and receive data between each pair of processors must be a contiguous
     *  part of the sendData or recvData.
     *
     *  The size of recvData must be recvPlan.totalQuantity().
     *  The size of sendData must be sendPlan.totalQuantity().
     */
    template<typename ValueType>
    tasking::SyncToken* exchangeByPlanAsync(
        ValueType recvVals[],
        const CommunicationPlan& recvPlan,
        const ValueType sendVals[],
        const CommunicationPlan& sendPlan ) const;

    /** Non-template version with coded ValueType is pure routine to be implemented by derived classes. */

    virtual tasking::SyncToken* exchangeByPlanAsyncImpl(
        void* recvData,
        const CommunicationPlan& recvPlan,
        const void* sendData,
        const CommunicationPlan& sendPlan,
        const common::ScalarType stype ) const = 0;

    /* @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myVals values that I receive
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allVals values for all processors (size must be n * size() )
     */
    template<typename ValueType>
    void scatter( ValueType myVals[], const IndexType n, const PartitionId root, const ValueType allVals[] ) const;

    /** Non-template version with coded ValueType is pure routine to be implemented by derived classes. */

    virtual void scatterImpl( void* myVals, const IndexType n, const PartitionId root,
                              const void* allVals, const common::ScalarType stype ) const = 0;

    /* @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myVals values that I receive
     *  @param[in]    n      number of elements in vector myVals
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    allVals values for all processors (size must be sum(sizes) )
     *  @param[in]    sizes   number of total values for all processors
     */
    template<typename ValueType>
    void scatterV( ValueType myVals[], const IndexType n, const PartitionId root, const ValueType allVals[], const IndexType sizes[] ) const;

    /* @brief Scatter of an array of values from root to all other processors.
     *
     *  @param[out]   myVals values that I receive
     *  @param[in]    root   processor that has values for all processors
     *  @param[in]    allVals values for all processors ( only root )
     *  @param[in]    sizes   number of values for each processor ( only root )
     */
    template<typename ValueType>
    void scatterVArray( hmemo::HArray<ValueType>& myVals, 
                        const PartitionId root, 
                        const hmemo::HArray<ValueType>& allVals, 
                        const hmemo::HArray<IndexType>& sizes ) const;

    virtual void scatterVImpl( void* myVals, const IndexType n, const PartitionId root,
                               const void* allVals, const IndexType sizes[], const common::ScalarType stype ) const = 0;

    /* @brief Gather of an array of values from all processors to root.
     *
     *  @param[out]   allVals values that I receive (size must be n * size() )
     *  @param[in]    n      number of elements in vector val
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myVals values that this processor contributes
     *  @param[in]    stype  specifies the data type used
     */
    template<typename ValueType>
    void gather( ValueType allVals[], const IndexType n, const PartitionId root, const ValueType myVals[] ) const;

    virtual void gatherImpl( void* allVals, const IndexType n, const PartitionId root, const void* myVals, const common::ScalarType stype ) const = 0;

    /* @brief Gather of an array of double values from all processors to root.
     *
     *  @param[out]   allVals values that I receive (size must be sum(sizes) )
     *  @param[in]    n      number of elements in myvals
     *  @param[in]    root   processor with values for all processors
     *  @param[in]    myVals values on this processor
     *  @param[in]    sizes  array with sizes from each processor
     *
     *  The array sizes must only be valid on the root processor. The
     *  output array allvals will only be valid on the root processor.
     */
    template<typename ValueType>
    void gatherV( ValueType allVals[], const IndexType n, const PartitionId root, const ValueType myVals[], const IndexType sizes[] ) const;

    virtual void gatherVImpl(
        void* allvals,
        const IndexType n,
        const PartitionId root,
        const void* myvals,
        const IndexType sizes[],
        const common::ScalarType stype ) const = 0;

    /* @brief Swap of an array with another processor.
     *
     * @param[in,out] val is the data array to be swapped
     * @param[in] n is the number of entries in array val
     * @param[in] partner is the rank of partition with which this partition swaps
     *
     * This method can also be used if partner is same as this processor.
     */
    template<typename ValueType>
    void swap( ValueType val[], const IndexType n, const PartitionId partner ) const;

    virtual void swapImpl(
        void* val,
        const IndexType n,
        PartitionId partner,
        const common::ScalarType stype ) const = 0;

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
     */

    template<typename ValueType>
    IndexType shift( ValueType newVals[], const IndexType newSize,
                     const ValueType oldVals[], const IndexType oldSize,
                     const int direction ) const;

    virtual IndexType shiftImpl(
        void* newVals,
        const IndexType newSize,
        const PartitionId source,
        const void* oldVals,
        const IndexType oldSize,
        const PartitionId dest,
        const common::ScalarType stype ) const = 0;

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
     */

    template<typename ValueType>
    tasking::SyncToken* shiftAsync(
        ValueType recvVals[],
        const ValueType sendVals[],
        const IndexType size,
        const int direction ) const;

    /** Pure method that provides the wanted functionality for each communicator */

    virtual tasking::SyncToken* shiftAsyncImpl(
        void* newVals,
        const PartitionId source,
        const void* oldVals,
        const PartitionId dest,
        const IndexType size,
        const common::ScalarType stype ) const = 0;

    /** All-to-all exchange of a single value between all processors.
     *
     * @tparam ValueType specifies the type of the elements for exchange
     *
     * @param[out] recvValues   will contain one value from each processor
     * @param[in]  sendValues   must contain one value for each processor
     *
     * recvValues and sendValues must both have a size of communicator size.
     * recvValues[i] on processor j contains sendValues[j] of processor i.
     */
    template<typename ValueType>
    void all2all( ValueType recvValues[], const ValueType sendValues[] ) const;

    /** Same routine as all2all but uses void pointers and codes the ValueType so it can
     *  become a virtual method that is provided by all derived communicator classes.
     */
    virtual void all2allImpl( 
        void* recvBuffer, 
        const void* sendBuffer,
        const common::ScalarType stype ) const = 0;

    /** All-to-all exchange of multiple values between all processors
     *
     *  @param[out] recvBuffer is an array where recvBuffer[i] is pointer to the data that will be received from processor i
     *  @param[in]  recvCount is an array where recvCount[i] is the number of values to receive from processor i
     *  @param[in]  sendBuffer is an array where sendBuffer[i] is pointer to the data to send to processor i
     *  @param[in]  sendCount is an array where sendCount[i] is the number of values to send to processor i
     *
     *  Here the values of recvCount must already be available. Before, it might have been computed by all2all( recvCount, sendCount ).
     */
    template<typename ValueType>
    void all2allv( ValueType* recvBuffer[], const IndexType recvCount[],
                   const ValueType* sendBuffer[], const IndexType sendCount[] ) const;

    /** Same routine but uses void pointers and codes the ValueType so it can
     *  become a virtual method that is provided by all derived communicator classes.
     */

    virtual void all2allvImpl(
        void* recvBuffer[], const IndexType recvCount[],
        const void* sendBuffer[], const IndexType sendCount[],
        const common::ScalarType stype ) const = 0;

    /** Pure method for bcast
     *
     *  @param[in,out] values  in on root, out on all other processors
     *  @param[in]     n    number of elements in vector val
     *  @param[in]     root processor with correct values of val
     *  @param[in]     stype codes the used data type of values
     */
    virtual void bcastImpl( void* values, const IndexType n, const PartitionId root, const common::ScalarType stype ) const = 0;

    /**************************************************************************************
     *                                                                                    *
     *  Other Virtual routines provided by derived class                                  *
     *                                                                                    *
     *************************************************************************************/

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

    /* @brief Broadcast a typed array from root to all other processors.
     *
     *  @tparam ValueType  is the data type
     *  @param[in,out] val  in on root, out on all other processors
     *  @param[in]     n    number of elements in vector val
     *  @param[in]     root processor with correct values of val
     *
     *  Note: implementation uses pure method bcastImpl without template argument
     */
    template<typename ValueType>
    void bcast( ValueType val[], const IndexType n, const PartitionId root ) const;

    /** Broadcast of a string.
     *
     *  This class provides one implementation, but derived classes might override it.
     */
    virtual void bcast( std::string&, const PartitionId root ) const;

    /* @brief Sum  up one single value from each partition to a global value.
     *
     *  @param[in] localValue  value on the calling partition
     *  @returns   global value, available for all partitions.
     *
     *  Implementation uses pure virtual method sumImpl (without template arguments)
     */
    template<typename ValueType>
    inline ValueType sum( const ValueType localValue ) const;

    template<typename ValueType>
    inline ValueType min( const ValueType localValue ) const;

    template<typename ValueType>
    inline ValueType max( const ValueType localValue ) const;

    /**  Sum values from all processes and distributes the result back to all processes.
     *
     *   @param[out] outValues array with the result values, same on all processes
     *   @param[in]  inValues individual contributions on each process
     *   @param[in]  n is the number of values in arrays inValues and outValues
     *   @param[in]  stype specifies the data type of the data
     */
    virtual void sumImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const = 0;

    /**  Find minimal values from all processes and distributes the result back to all processes. */

    virtual void minImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const = 0;

    /**  Find maximal values from all processes and distributes the result back to all processes. */

    virtual void maxImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const = 0;

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

    template<typename ValueType>
    void maxloc( ValueType& val, IndexType& location, const PartitionId root ) const;

    template<typename ValueType>
    void minloc( ValueType& val, IndexType& location, const PartitionId root ) const;

    virtual void maxlocImpl( void* val, IndexType* location, PartitionId root, const common::ScalarType stype ) const = 0;

    virtual void minlocImpl( void* val, IndexType* location, PartitionId root, const common::ScalarType stype ) const = 0;

    /** Scan values among the processor belonging to this communicator. */

    template<typename ValueType>
    ValueType scan( const ValueType localValue ) const;

    /** Inclusive scan, similiar to sum but each processor has partial sums */

    virtual void scanImpl( void* outValues, const void* inValues, const IndexType n, const common::ScalarType stype ) const = 0;

    /** Default implementation */

    template<typename ValueType>
    ValueType scanDefault( const ValueType localValue ) const;

    /**
     *  Predicate that returns true if (derived) Communicator class supports minloc/maxlocImpl for a
     *  given value type.
     *
     *  @param[in] vType specifies the type of the reduction array values
     *  @param[in] iType specifies the index type used
     */
    virtual bool supportsLocReduction( const common::ScalarType vType, const common::ScalarType iType ) const = 0;

    /** Default implementation for maxloc that uses a gather operation instead of reduction. */

    template<typename ValueType>
    void maxlocDefault( ValueType& val, IndexType& location, const PartitionId root ) const;

    /** Default implementation for minloc that uses a gather operation instead of reduction. */

    template<typename ValueType>
    void minlocDefault( ValueType& val, IndexType& location, const PartitionId root ) const;

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

    void computeOwners(
        hmemo::HArray<PartitionId>& owners,
        const Distribution& distribution,
        const hmemo::HArray<IndexType>& requiredIndexes ) const;

    /**************************************************************************************
     *                                                                                    *
     *  Communication routines for HArrays for convenience                             *
     *                                                                                    *
     *************************************************************************************/

    /** Broadcast of a heterogeneous array */

    template<typename ValueType>
    void bcastArray( hmemo::HArray<ValueType>& array, const PartitionId root ) const;

    /** Broadcast of a heterogeneous array with known size
     *
     *  If the size of the array is known, an additional broadcast of the size is not required.
     */

    template<typename ValueType>
    void bcastArray( hmemo::HArray<ValueType>& array, const IndexType n, const PartitionId root ) const;

    /** exchangeByPlan for HArrays instead of usual array */

    template<typename ValueType>
    void exchangeByPlan(
        hmemo::HArray<ValueType>& recvArray,
        const CommunicationPlan& recvPlan,
        const hmemo::HArray<ValueType>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** exchangeByPlan as function with result argument for convenience. */

    template<typename ValueType>
    hmemo::HArray<ValueType> exchangeByPlanF(
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

    template<typename ValueType>
    void joinArray( hmemo::HArray<ValueType>& global, const hmemo::HArray<ValueType>& local ) const;

    /** @brief Asychronous shift on LAMA arrays.
     *
     *  @tparam     ValueType   stands for the data type of the arrays
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

    /** @brief Sum up values of an array
     *
     *  @tparam        ValueType   specifies the type of the array data
     *  @param[in,out] array       in are the local values, out are the global summed values
     *
     *  Note: All partitions must have the same size for the array
     */
    template<typename ValueType>
    void sumArray( hmemo::HArray<ValueType>& array ) const;

    /** Override routine of base class Printable. */

    virtual void writeAt( std::ostream& stream ) const;

    /** @brief Getter for the type of a communicator. */

    inline CommunicatorType getType() const;

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

    /**
     *  @brief This routine provides a collective / concurrent file that allows parallel I/O
     */
    virtual std::unique_ptr<class CollectiveFile> collectiveFile() const = 0;

    /** Read in the environment variable SCAI_NP for user processor array.
     *
     * @param[out] userProcArray specifies the user processor array.
     *
     */
    static void getUserProcArray( PartitionId userProcArray[3] );

    void setSeed( int seed ) const;

    /**
     *  @brief  Function that returns the name of node on which this communicator runs
     */
    const char* getNodeName() const;

    /**
     *  @brief Function that returns a unique id of the node on which this process is runing.
     *
     *  The unique id cal also be considered as the rank of this node among all other nodes.
     *  The rank of a node is given by sorting by the ranks of those processors with node rank equal 0.
     */
    int getNodeId() const;

protected:

    /** Constructor of abstract classes are always protected. */

    Communicator( const CommunicatorType& type );

    /** This method determines node rank and node size by comparing the names. */

    void setNodeData();

    void setSizeAndRank( PartitionId size, PartitionId rank );

    /** Get the processor name.
     *
     *  @param[out] name is the processor name, allocated with at least maxProcessorName
     */
    virtual void getProcessorName( char* name ) const = 0;

    /** Pure method that returns the maximal length of string set by getProcessorName */

    virtual size_t maxProcessorName() const = 0;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Shift implementation for direction == 0, just copies values. */

    template<typename ValueType>
    IndexType shift0( ValueType newVals[], const IndexType newSize,
                      const ValueType oldVals[], const IndexType oldSize ) const;

private:

    CommunicatorType mCommunicatorType; //!< type of this communicator

    PartitionId mRank; //!< rank of this processor

    PartitionId mSize; //!< number of processors in this communicator

    PartitionId mNodeRank; //!< rank of this processor on its node

    PartitionId mNodeSize; //!< number of processors on same node

    PartitionId mIdNode;  //!< unique id of this node

    PartitionId mNumNodes;  //!< total number of nodes used

    std::unique_ptr<char[]> mNodeName;

    mutable int mSeed;
};

/* -------------------------------------------------------------------------- */
/*  Implementation of inline methods                                          */
/* -------------------------------------------------------------------------- */

PartitionId Communicator::getSize() const
{
    return mSize;
}

/* -------------------------------------------------------------------------- */

PartitionId Communicator::getRank() const
{
    return mRank;
}

/* -------------------------------------------------------------------------- */

PartitionId Communicator::getNodeSize() const
{
    return mNodeSize;
}

/* -------------------------------------------------------------------------- */

PartitionId Communicator::getNodeRank() const
{
    return mNodeRank;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Communicator::sum( ValueType localValue ) const
{
    ValueType globalValue;

    // general routine uses typeless pointers void*, type is coded via ScalarType

    sumImpl( &globalValue, &localValue, 1, common::TypeTraits<ValueType>::stype );

    return globalValue;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Communicator::min( ValueType localValue ) const
{
    ValueType globalValue;

    // general routine uses typeless pointers void*, type is coded via ScalarType

    minImpl( &globalValue, &localValue, 1, common::TypeTraits<ValueType>::stype );

    return globalValue;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Communicator::max( ValueType localValue ) const
{
    ValueType globalValue;

    // general routine uses typeless pointers void*, type is coded via ScalarType

    maxImpl( &globalValue, &localValue, 1, common::TypeTraits<ValueType>::stype );

    return globalValue;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::bcast( ValueType val[], const IndexType n, const PartitionId root ) const
{
    SCAI_ASSERT_LT_ERROR( root, getSize(), *this << ": Illegal root " << root )

    // general pure routine uses typeless pointers void*, type is coded via ScalarType

    bcastImpl( val, n, root, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::gather( ValueType allVals[], const IndexType n, const PartitionId root, const ValueType myVals[] ) const
{
    gatherImpl( allVals, n, root, myVals, common::TypeTraits<ValueType>::stype );
}

template<typename ValueType>
void Communicator::scatter( ValueType myVals[], const IndexType n, const PartitionId root, const ValueType allVals[] ) const
{
    scatterImpl( myVals, n, root, allVals, common::TypeTraits<ValueType>::stype );
}

template<typename ValueType>
void Communicator::gatherV( ValueType allVals[], const IndexType n, const PartitionId root, const ValueType myVals[], const IndexType sizes[] ) const
{
    gatherVImpl( allVals, n, root, myVals, sizes, common::TypeTraits<ValueType>::stype );
}

template<typename ValueType>
void Communicator::scatterV( ValueType myVals[], const IndexType n, const PartitionId root, const ValueType allVals[], const IndexType sizes[] ) const
{
    scatterVImpl( myVals, n, root, allVals, sizes, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* Communicator::exchangeByPlanAsync(
    ValueType recvVals[],
    const CommunicationPlan& recvPlan,
    const ValueType sendVals[],
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvVals, recvPlan, sendVals, sendPlan, common::TypeTraits<ValueType>::stype );
}

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
void Communicator::swap( ValueType val[], const IndexType n, const PartitionId partner ) const
{
    swapImpl( val, n, partner, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void Communicator::exchangeByPlan(
    ValueType recvVals[],
    const CommunicationPlan& recvPlan,
    const ValueType sendVals[],
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvVals, recvPlan, sendVals, sendPlan, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType Communicator::shift(
    ValueType newVals[], const IndexType newSize,
    const ValueType oldVals[], const IndexType oldSize,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        // source and dest are the same processor

        return shift0( newVals, newSize, oldVals, oldSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );

    return shiftImpl( newVals, newSize, source, oldVals, oldSize, dest, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* Communicator::shiftAsync(
    ValueType recvVals[],
    const ValueType sendVals[],
    const IndexType size,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        // source and dest are the same processor

        shift0( recvVals, size, sendVals, size );
        return new tasking::NoSyncToken();
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );

    return shiftAsyncImpl( recvVals, source, sendVals, dest, size, common::TypeTraits<ValueType>::stype );
}

/* -------------------------------------------------------------------------- */

CommunicatorType Communicator::getType() const
{
    return mCommunicatorType;
}

} /* end namespace dmemo */

} /* end namespace scai */

/* -------------------------------------------------------------------------- */

/** Macro that defines a new current communicator for the actual scope
 *
 *  \code
 *      auto comm = Communicator::getCommunicatorPtr();
 *      PartitionId color = ...;
 *      auto commTask = comm->split( color );
 *      {
 *          SCAI_DMEMO_TASK( commTask )
 *          auto dist = std::make_shared<BlockDistribution>( n );  // distributes onto commTask
 *          ...
 *      }
 *  \endcode
 */
#define SCAI_DMEMO_TASK( comm ) scai::dmemo::ScopedCommunicatorRecord SCAI_Comm_( comm );

