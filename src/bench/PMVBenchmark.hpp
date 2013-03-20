/**
 * @file PMVBenchmark.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Benchmark for parallel matrix-vector multiplication with different settings.
 * @author Thomas Brandes, Jiri Kraus
 * @date 11.05.2011
 * $Id$
 */
#ifndef LAMA_PMV_BENCHMARK_HPP_
#define LAMA_PMV_BENCHMARK_HPP_

#include <bench/LAMAMPIBenchmark.hpp>
#include <bench/LAMAInputSet.hpp>
#include <bench/LAMAInputSetComplexityVisitor.hpp>

#include <lama/DenseVector.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>
#include <lama/tracing.hpp>

#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/GenBlockDistribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <lama/norm/MaxNorm.hpp>

#ifdef LAMA_BUILD_CUDA
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <map>
#include <string>

using lama::LAMAInterfaceRegistry;
using lama::Context;
using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PMVBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PMVBenchmark();

    PMVBenchmark( const std::string& arguments );

    PMVBenchmark( const PMVBenchmark<MatrixType>& other );

    virtual ~PMVBenchmark();

    virtual std::auto_ptr<bf::Benchmark> copy() const;

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    virtual const std::string& getId() const;

protected:
    virtual void initialize();
    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

    using bf::Benchmark::mName;
    using bf::Benchmark::mGId;

    using LAMAMPIBenchmark::mComm;

private:

    static const LAMAInputSetComplexityVisitor::Group& group();

    static const std::string& sid();

    std::auto_ptr<MatrixType> mMatrixA;
    std::auto_ptr<lama::DenseVector<ValueType> > mVectorX;
    std::auto_ptr<lama::DenseVector<ValueType> > mVectorY;

    lama::DistributionPtr mDistribution;
    lama::ContextPtr mLocalContext;
    lama::ContextPtr mHaloContext;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
};

template<typename MatrixType>
const std::string& PMVBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PMVBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mLocalContext(
        lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mLocalContext(
        lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark( const PMVBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mLocalContext( lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
        lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
            0 ), mNumProcessedBytesDouble( 0 )
{

}

template<typename MatrixType>
PMVBenchmark<MatrixType>::~PMVBenchmark()
{
    // Note: mMatrixA, mVectorX, mVectorY will be freed
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PMVBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PMVBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PMVBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PMVBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PMVBenchmark<MatrixType>::getId() const
{
    return id();
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::initialize()
{
    //device initialization + comunicator

    map<string,string> tokens;

    getConfig( tokens );

    int noThreads = 1;
    int devNo = -1; // default device

    if( tokens.count( "CUDA" ) > 0 )
    {
        istringstream tokenStream( tokens["CUDA"] );
        tokenStream >> devNo;
    }

    std::ostringstream idStream;
    idStream << " (Proc " << mComm->getRank() << " LOCAL=";
    if( tokens["LOCAL"] == "CUDA" )
    {
        mLocalContext = lama::ContextFactory::getContext( Context::CUDA, devNo );
        idStream << "CUDA:";
    }
    else
    {
        mLocalContext = lama::ContextFactory::getContext( Context::Host );
        idStream << "HOST:";
    }
    idStream << "HALO=";
    if( tokens["HALO"] == "CUDA" )
    {
        mHaloContext = lama::ContextFactory::getContext( Context::CUDA, devNo );
        idStream << "CUDA:";
    }
    else if( tokens["HALO"] == "HOST" )
    {
        mHaloContext = lama::ContextFactory::getContext( Context::Host );
        idStream << "HOST:";
    }
    else
    {
        mHaloContext = mLocalContext;
        idStream << "LOCAL:";
    }
    if( tokens.count( "THREADS" ) > 0 )
    {
        istringstream tokenStream( tokens["THREADS"] );
        tokenStream >> noThreads;
    }
    else
    {
        noThreads = 1;
    }
    idStream << "THREADS=" << noThreads << ":";
    if( mLocalContext->getType() == lama::Context::CUDA || mHaloContext->getType() == lama::Context::CUDA )
    {
        idStream << "DEVICE=" << devNo << ":";
    }
    idStream << "COMM=";
    if( tokens["COMM"] == "SYNC" )
    {
        mCommunicationKind = lama::Matrix::SYNCHRONOUS;
        idStream << "SYNC:";
    }
    else
    {
        mCommunicationKind = lama::Matrix::ASYNCHRONOUS;
        idStream << "ASYNC:";
    }
    float weight = 1.0;
    if( tokens.count( "W" ) > 0 )
    {
        istringstream tokenStream( tokens["W"] );
        tokenStream >> weight;
    }
    idStream << "W=" << weight;

    //mName += idStream.str();

    //TODO: Aggregate Name at Proc 0 to print it in output
    printf( " Name: %s )\n", idStream.str().c_str() );
    omp_set_num_threads( noThreads );
#ifdef LAMA_BUILD_CUDA
    if( LAMAInterfaceRegistry::getRegistry().hasInterface( lama::Context::CUDA ) )
    {
        // CUDA context is supported, enable CUDAHost for Host if CUDA is used

        if( mLocalContext->getType() == lama::Context::CUDA )
        {
            // support fast memory transfer Host->CUDA
            lama::CUDAHostContextManager::setAsCurrent( mLocalContext );
        }
        else if( mHaloContext->getType() == lama::Context::CUDA )
        {
            // support fast memory transfer Host->CUDA
            lama::CUDAHostContextManager::setAsCurrent( mHaloContext );
        }

        // Note: if local and halo context use both but different CUDA devices, fast memory
        //       transfer is only supported to the device for local computations
    }
#endif //LAMA_BUILD_CUDA
    LAMA_LOG_INFO( logger, "get input set " << mInputSetId );

    const LAMAInputSet& inputSet = bf::InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    const lama::DenseVector<double>& inputX = inputSet.getX();
    const lama::CSRSparseMatrix<double>& inputA = inputSet.getA();

    LAMA_LOG_INFO( logger, "input matrix A = " << inputA );

    // mDistribution = lama::DistributionPtr ( new lama::GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );

    mDistribution = inputA.getDistributionPtr();

    LAMA_LOG_INFO( logger,
                   "General block distribution of matrix with weight " << weight << ", local size = " << mDistribution->getLocalSize() );

    mVectorX.reset( new lama::DenseVector<ValueType>( inputX, mDistribution ) );
    mVectorY.reset( new lama::DenseVector<ValueType>( mDistribution, 0.0 ) );

    LAMA_LOG_DEBUG( logger, "convert and redistribute the input matrix " << inputA );

    mMatrixA.reset( new MatrixType( inputA ) );

    mMatrixA->redistribute( mDistribution, mDistribution );

    LAMA_LOG_DEBUG( logger, "matrix " << *mMatrixA << " now available for MV" );

    mMatrixA->setContext( mLocalContext, mHaloContext );
    mMatrixA->setCommunicationKind( mCommunicationKind );
    mMatrixA->setContext( mLocalContext, mHaloContext );

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << mLocalContext << ", halo at " << mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::setUp()
{
    mVectorX->prefetch( mLocalContext );
    mMatrixA->prefetch();
    mVectorX->wait();
    mMatrixA->wait();

    LAMA_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : X, A locat at " << mLocalContext << " and A halo at " << mHaloContext );
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::execute()
{
    LAMA_REGION( "execute" );
    lama::Vector& result = *mVectorY;
    const lama::Matrix& A = *mMatrixA;
    const lama::Vector& x = *mVectorX;
    result = A * x;
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::tearDown()
{
    lama::ContextPtr host = lama::ContextFactory::getContext( Context::Host );
    mVectorY->prefetch( host );
    mVectorY->wait();
    LAMA_LOG_INFO( logger, "tearDown done, Y at Host" );
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::shutdown()
{
    const LAMAInputSet& inputSet = bf::InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    //TODO: Complexity Calculation needs to be ported
    if( ( typeid(MatrixType) == typeid(lama::DenseMatrix<float>) )
            || ( typeid(MatrixType) == typeid(lama::DenseMatrix<double>) ) )
    {
        LAMAInputSetComplexityVisitor::getMVComplexity( inputSet.getDenseA(), mNumFloatingPointOperations,
                mNumProcessedBytesFloat, mNumProcessedBytesDouble );
    }
    else
    {
        LAMAInputSetComplexityVisitor::getMVComplexity( inputSet.getA(), mNumFloatingPointOperations,
                mNumProcessedBytesFloat, mNumProcessedBytesDouble );
    }

    const lama::DenseVector<double>& result = inputSet.getY();

    // set the maximal difference that is allowed, depends on ValueType

    ValueType maxDiff = 1E-4f;

    if( typeid(ValueType) == typeid(double) )
    {
        maxDiff = maxDiff * maxDiff; // double precision
    }

    try
    {
        const lama::DenseVector<ValueType>& computedResult = *mVectorY;

        LAMA_LOG_INFO( logger, "computed result: " << computedResult );

        // Note: there might be type conversion from double to float

        lama::DenseVector<ValueType> correctResult( result );

        LAMA_LOG_INFO( logger, "result ( by inputSet ): " << result );
        LAMA_LOG_INFO( logger, "correct result: " << correctResult );

#if defined( LAMA_LOG_TRACE_ENABLED)
        for ( int i = 0; i < correctResult.size(); i++)
        {
            ValueType inputValue = lama::cast<ValueType>( ( *mVectorX )( i ) );
            ValueType correctValue = lama::cast<ValueType>( correctResult(i) );
            ValueType computedValue = lama::cast<ValueType>( computedResult(i) );
            LAMA_LOG_TRACE( logger, i << ": correct = " << correctValue
                            << ", computed = " << computedValue
                            << ", input = " << inputValue );
        }
#endif

        lama::Vector& diff = correctResult;

        diff = computedResult - correctResult;

        lama::Scalar diffNorm = lama::maxNorm( diff );

        LAMA_LOG_INFO( logger, "max diff = " << diffNorm );

        LAMA_CHECK_BENCHMARK( diffNorm.getValue<ValueType>() < maxDiff );

    }
    catch( bf::BFException& e )
    {
        std::stringstream message;
        message << e.what() << " std::fabs( result[i] - computedResultValue ) is bigger than " << maxDiff << std::endl;
        throw bf::BFException( message.str() );
    }

    mMatrixA.reset();
    mVectorX.reset();
    mVectorY.reset();
}

template<typename MatrixType>
CounterType PMVBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return mNumFloatingPointOperations;
}

template<typename MatrixType>
CounterType PMVBenchmark<MatrixType>::getProcessedBytes() const
{
    typedef typename MatrixType::ValueType ValueType;

    if( sizeof(ValueType) == sizeof(float) )
    {
        return mNumProcessedBytesFloat;
    }
    else if( sizeof(ValueType) == sizeof(double) )
    {
        return mNumProcessedBytesDouble;
    }
    return 0;
}

#endif // LAMA_PMV_BENCHMARK_HPP_
