/**
 * @file NoCommunicator.hpp
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
 * @brief Implementation of a communicator class for non-distributed (replicated) objects.
 * @author Thomas Brandes
 * @date 15.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_NO_COMMUNICATOR_HPP_
#define LAMA_NO_COMMUNICATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/CRTPCommunicator.hpp>

namespace lama
{

/** The class NoCommunicator stands for objects that are replicated on each
 *  partition or processor.
 */

class LAMA_DLL_IMPORTEXPORT NoCommunicator: public CRTPCommunicator<NoCommunicator>
{

    friend class CRTPCommunicator<NoCommunicator> ;

public:

    NoCommunicator();

    virtual ~NoCommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    virtual PartitionId getSize() const;

    virtual PartitionId getRank() const;

    virtual PartitionId getNodeSize() const;

    virtual PartitionId getNodeRank() const;

    virtual void all2all( IndexType recvValues[], const IndexType sendValues[] ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private    :

    // Implementation methods are all private, but CRTPCommunicator is a friend class

    template<typename T>
    IndexType shiftImpl(
                    T newvals[],
                    const IndexType newSize,
                    const PartitionId source,
                    const T oldVals[],
                    const IndexType oldSize,
                    const PartitionId dest ) const;

    template<typename T>
    SyncToken* shiftAsyncImpl(
                    T newvals[],
                    const PartitionId source,
                    const T oldVals[],
                    const PartitionId dest,
                    const IndexType size ) const;

    template<typename T>
    void bcastImpl( T val[], const IndexType n, const PartitionId root ) const;

    template<typename T>
    void scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const;

    template<typename T>
    void scatterVImpl(
                    T myvals[],
                    const IndexType n,
                    const PartitionId root,
                    const T allvals[],
                    const IndexType sizes[] ) const;

    template<typename T>
    void gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const;

    template<typename T>
    void gatherVImpl(
                    T allvals[],
                    const IndexType n,
                    const PartitionId root,
                    const T myvals[],
                    const IndexType sizes[] ) const;

    template<typename T>
    T sumImpl( const T value ) const;

    template<typename T>
    T maxImpl( const T value ) const;

    template<typename T>
    T minImpl( const T value ) const;

    template<typename T>
    void maxlocImpl( T& val, IndexType& location, const PartitionId root ) const;

    template<typename T>
    void swapImpl( T val[], const IndexType n, const PartitionId partner ) const;

    // common implementation for self exchange, uses size of datatype

    template<typename T>
    void exchangeByPlanImpl(
                    T recvData[],
                    const CommunicationPlan& recvPlan,
                    const T sendData[],
                    const CommunicationPlan& sendPlan ) const;

    template<typename T>
    SyncToken* exchangeByPlanAsyncImpl(
                    T recvData[],
                    const CommunicationPlan& recvPlan,
                    const T sendData[],
                    const CommunicationPlan& sendPlan ) const;

    virtual ContextPtr getCommunicationContext( const _LAMAArray& array ) const;
};

}

#endif // LAMA_NO_COMMUNICATOR_HPP_
