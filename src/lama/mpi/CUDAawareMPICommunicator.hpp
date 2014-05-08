/**
 * @file CUDAawareMPICommunicator.hpp
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
 * @brief CUDAawareMPICommunicator.hpp
 * @author Lauretta Schubert
 * @date 11.03.2014
 * @since 1.0.1
 */

#ifndef LAMA_MPI_CUDAAWARECOMMUNICATOR_HPP_
#define LAMA_MPI_CUDAAWARECOMMUNICATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/mpi/MPICommunicator.hpp>

namespace lama
{

/** Communicator class that implements communication and data exchange via MPI with CUDAawareness.
 *
 *  MPI_Init is called in the constructor, MPI_Finalize is called in the destructor.
 */

class LAMA_DLL_IMPORTEXPORT CUDAawareMPICommunicator: public MPICommunicator
{

// Only MPICommunicatorManager is allowed to create MPI communicator

    friend class CUDAawareMPICommunicatorManager;

public:
    virtual ~CUDAawareMPICommunicator();

    /** Exchange of LAMAArrays. */

    template<typename T>
    void exchangeByPlan(
        LAMAArray<T>& recvArray,
        const CommunicationPlan& recvPlan,
        const LAMAArray<T>& sendArray,
        const CommunicationPlan& sendPlan ) const;

    /** Asynchronous exchange of LAMAArrays. */

    template<typename T>
    SyncToken* exchangeByPlanAsync(
        LAMAArray<T>& recvArray,
        const CommunicationPlan& recvPlan,
        const LAMAArray<T>& sendArray,
        const CommunicationPlan& sendPlan ) const;

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

    virtual void writeAt( std::ostream& stream ) const;

private:

    CUDAawareMPICommunicator( int& argc, char** & argv );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static const int defaultTag;
};

template<typename T>
void CUDAawareMPICommunicator::exchangeByPlan(
    LAMAArray<T>& recvArray,
    const CommunicationPlan& recvPlan,
    const LAMAArray<T>& sendArray,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_ERROR(
        sendArray.size() == sendPlan.totalQuantity(),
        "Send array has size " << sendArray.size() << ", but send plan requires " << sendPlan.totalQuantity() << " entries" )

    IndexType recvSize = recvPlan.totalQuantity();

    WriteAccess<T> recvData( recvArray, recvArray.getValidContext() );
    ReadAccess<T> sendData( sendArray, sendArray.getValidContext() );

    recvData.clear();
    recvData.resize( recvSize );

    exchangeByPlan( recvData.get(), recvPlan, sendData.get(), sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename T>
SyncToken* CUDAawareMPICommunicator::exchangeByPlanAsync(
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

    boost::shared_ptr<WriteAccess<T> > recvData( new WriteAccess<T>( recvArray, recvArray.getValidContext() ) );
    boost::shared_ptr<ReadAccess<T> > sendData( new ReadAccess<T>( sendArray, sendArray.getValidContext() ) );

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

#endif // LAMA_MPI_CUDAAWARECOMMUNICATOR_HPP_
