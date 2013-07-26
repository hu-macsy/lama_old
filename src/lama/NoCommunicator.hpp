/**
 * @file NoCommunicator.hpp
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
#include <lama/Communicator.hpp>

// boost
#include <boost/weak_ptr.hpp>

namespace lama
{

/** The class NoCommunicator stands for objects that are replicated on each
 *  partition or processor.
 */

class LAMA_DLL_IMPORTEXPORT NoCommunicator: public Communicator
{
public:

    NoCommunicator();

    virtual ~NoCommunicator();

    virtual bool isEqual( const Communicator& other ) const;

    virtual ThreadSafetyLevel getThreadSafetyLevel() const;

    virtual PartitionId getSize() const;

    virtual PartitionId getRank() const;

    virtual PartitionId getNodeSize() const;

    virtual PartitionId getNodeRank() const;

    virtual void all2all( int* recvValues, const int* sendValues ) const;

    virtual void exchangeByPlan(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual void exchangeByPlan(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual void exchangeByPlan(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual SyncToken* exchangeByPlanAsync(
        int* const recvData,
        const CommunicationPlan& recvPlan,
        const int* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual SyncToken* exchangeByPlanAsync(
        float* const recvData,
        const CommunicationPlan& recvPlan,
        const float* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual SyncToken* exchangeByPlanAsync(
        double* const recvData,
        const CommunicationPlan& recvPlan,
        const double* const sendData,
        const CommunicationPlan& sendPlan ) const;

    virtual IndexType shiftImpl(
        double targetVals[],
        const IndexType targetSize,
        const double sourceVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual IndexType shiftImpl(
        float targetVals[],
        const IndexType targetSize,
        const float sourceVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual IndexType shiftImpl(
        int targetVals[],
        const IndexType targetSize,
        const int sourceVals[],
        const IndexType oldSize,
        const int direction ) const;

    virtual void bcast( double val[], const IndexType n, const PartitionId root ) const;
    virtual void bcast( float val[], const IndexType n, const PartitionId root ) const;
    virtual void bcast( int val[], const IndexType n, const PartitionId root ) const;
    virtual void bcast( std::string& val, const PartitionId root ) const;

    virtual void scatter( double myvals[], const IndexType n, const PartitionId root, const double allvals[] ) const;
    virtual void scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const;
    virtual void scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const;

    virtual void scatter(
        double myvals[],
        const IndexType n,
        const PartitionId root,
        const double allvals[],
        const IndexType sizes[] ) const;
    virtual void scatter(
        float myvals[],
        const IndexType n,
        const PartitionId root,
        const float allvals[],
        const IndexType sizes[] ) const;
    virtual void scatter(
        int myvals[],
        const IndexType n,
        const PartitionId root,
        const int allvals[],
        const IndexType sizes[] ) const;

    virtual void gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const;
    virtual void gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const;
    virtual void gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const;

    virtual void gather(
        double allvals[],
        const IndexType n,
        const PartitionId root,
        const double myvals[],
        const IndexType sizes[] ) const;
    virtual void gather(
        float allvals[],
        const IndexType n,
        const PartitionId root,
        const float myvals[],
        const IndexType sizes[] ) const;
    virtual void gather(
        int allvals[],
        const IndexType n,
        const PartitionId root,
        const int myvals[],
        const IndexType sizes[] ) const;

    virtual float sum( const float value ) const;

    virtual double sum( const double value ) const;

    virtual int sum( const int value ) const;

    virtual size_t sum( const size_t value ) const;

    virtual float min( const float value ) const;

    virtual float max( const float value ) const;

    virtual double min( const double value ) const;

    virtual double max( const double value ) const;

    virtual int min( const int value ) const;

    virtual int max( const int value ) const;

    virtual void maxloc( double& val, int& location, const PartitionId root ) const;
    virtual void maxloc( float& val, int& location, const PartitionId root ) const;
    virtual void maxloc( int& val, int& location, const PartitionId root ) const;

    virtual void swap( double val[], const IndexType n, const PartitionId partner ) const;
    virtual void swap( float val[], const IndexType n, const PartitionId partner ) const;
    virtual void swap( int val[], const IndexType n, const PartitionId partner ) const;

    virtual void gather( std::vector<IndexType>& values, IndexType value ) const;

    virtual void gather( std::vector<float>& values, float value ) const;

    virtual void synchronize() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    // common implementation for self exchange, uses size of datatype

    void exchangeByPlanImpl(
        void* const recvData,
        const CommunicationPlan& recvPlan,
        const void* const sendData,
        const CommunicationPlan& sendPlan,
        int elemSize ) const;

    virtual ContextPtr getCommunicationContext() const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

#endif // LAMA_NO_COMMUNICATOR_HPP_
