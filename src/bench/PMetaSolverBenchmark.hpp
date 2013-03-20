/**
 * @file PMetaSolverBenchmark.hpp
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
 * @brief Benchmark for parallel MetaSolver implementation with different settings.
 * @author Jiri Kraus
 * @date 25.08.2011
 * $Id$
 */
#ifndef LAMA_LAMAMETASOLVERBENCHMARK_HPP_
#define LAMA_LAMAMETASOLVERBENCHMARK_HPP_

#include <sstream>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>

#include <lama/matrix/DenseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>

#include <lama/distribution/GenBlockDistribution.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <lama/solver/MetaSolver.hpp>

#ifdef LAMA_BUILD_CUDA
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <bench/LAMAMPIBenchmark.hpp>
#include <bench/LAMAInputSet.hpp>
#include <bench/LAMAInputSetComplexityVisitor.hpp>

using lama::LAMAInterfaceRegistry;
using lama::Context;
using lama::MetaSolver;
using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PMetaSolverBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PMetaSolverBenchmark();

    PMetaSolverBenchmark( const std::string& arguments );

    PMetaSolverBenchmark( const PMetaSolverBenchmark<MatrixType>& other );

    virtual ~PMetaSolverBenchmark();

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

    MatrixType* mMatrixA;

    lama::DenseVector<ValueType>* mSolution;
    lama::DenseVector<ValueType>* mRhs;

    lama::MetaSolver* mMetaSolver;

    std::string mFileName;

    lama::DistributionPtr mDistribution;
    lama::ContextPtr mLocalContext;
    lama::ContextPtr mHaloContext;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
};

template<typename MatrixType>
const std::string& PMetaSolverBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PMetaSolverBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PMetaSolverBenchmark<MatrixType>::PMetaSolverBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mMatrixA( 0 ), mSolution(
        0 ), mRhs( 0 ), mLocalContext( lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMetaSolverBenchmark<MatrixType>::PMetaSolverBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mMatrixA(
        0 ), mSolution( 0 ), mRhs( 0 ), mLocalContext(
            lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                    0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMetaSolverBenchmark<MatrixType>::PMetaSolverBenchmark( const PMetaSolverBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mMatrixA( 0 ), mSolution( 0 ), mRhs( 0 ), mLocalContext(
        lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{

}

template<typename MatrixType>
PMetaSolverBenchmark<MatrixType>::~PMetaSolverBenchmark()
{
    if( mMatrixA != 0 )
    {
        delete mMatrixA;
    }
    mMatrixA = 0;
    if( mSolution != 0 )
    {
        delete mSolution;
    }
    mSolution = 0;
    if( mRhs != 0 )
    {
        delete mRhs;
    }
    mRhs = 0;
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PMetaSolverBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PMetaSolverBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PMetaSolverBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PMetaSolverBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PMetaSolverBenchmark<MatrixType>::getId() const
{
    return id();
}

template<typename MatrixType>
void PMetaSolverBenchmark<MatrixType>::initialize()
{
    LAMA_LOG_INFO( logger, "initialize" );

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
        if( mLocalContext->getType() == Context::CUDA )
        {
            mHaloContext = mLocalContext;
        }
        else
        {
            mHaloContext = lama::ContextFactory::getContext( Context::CUDA, devNo );
        }
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

    LAMA_LOG_INFO( logger, "Set #threads = " << noThreads << ", omp_get_max_threads() = " << omp_get_max_threads() );

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

    LAMA_LOG_INFO( logger, "get matrix from input set" );

    const lama::CSRSparseMatrix<double>& inputA = inputSet.getA();

    LAMA_LOG_INFO( logger, "input matrix A = " << inputA );

    // mDistribution = lama::DistributionPtr ( new lama::GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );

    mDistribution = inputA.getDistributionPtr();

    // redistribute input matrix, is local copy in case of same distribution

    mMatrixA = new MatrixType( inputA );

    mMatrixA->redistribute( mDistribution, mDistribution );

    mMatrixA->setCommunicationKind( mCommunicationKind );
    mMatrixA->setContext( mLocalContext, mHaloContext );

    // redistribute vector Y of input set

    mRhs = new lama::DenseVector<ValueType>( inputSet.getY(), mDistribution );

    mSolution = new lama::DenseVector<ValueType>( mDistribution, 0.0 );

    mRhs->setContext( mLocalContext );
    mSolution->setContext( mLocalContext );

    /*
     * TODO  non static solution
     */
    mFileName = "config.txt";

    mMetaSolver = new MetaSolver( "MetaSolver solver", mFileName );

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << mLocalContext << ", halo at " << mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PMetaSolverBenchmark<MatrixType>::setUp()
{
    mMetaSolver->initialize( *mMatrixA );

    mSolution->prefetch( mLocalContext );
    mRhs->prefetch( mLocalContext );
    mMatrixA->prefetch(); // prefetch local storage to local and halo storage to halo context
    mSolution->wait();
    mRhs->wait();
    mMatrixA->wait();

    LAMA_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : Solution, Rhs and  A local at " << mLocalContext << " and A halo at " << mHaloContext );
}

template<typename MatrixType>
void PMetaSolverBenchmark<MatrixType>::execute()
{
    LAMA_LOG_INFO( logger, "execute benchmark" );

    mMetaSolver->solve( *mSolution, *mRhs );
}

template<typename MatrixType>
void PMetaSolverBenchmark<MatrixType>::tearDown()
{
    LAMA_LOG_INFO( logger, "tear down" );

    lama::ContextPtr host = lama::ContextFactory::getContext( Context::Host );

    mSolution->prefetch( host );
    mSolution->wait();

    LAMA_LOG_INFO( logger, "tearDown done, Y at Host" );
}

template<typename MatrixType>
void PMetaSolverBenchmark<MatrixType>::shutdown()
{
    LAMA_LOG_INFO( logger, "shutdown benchmark" );

    if( mMetaSolver != 0 )
    {
        delete mMetaSolver;
    }

    if( mMatrixA != 0 )
    {
        delete mMatrixA;
    }
    mMatrixA = 0;
    if( mSolution != 0 )
    {
        delete mSolution;
    }
    mSolution = 0;
    if( mRhs != 0 )
    {
        delete mRhs;
    }
    mRhs = 0;

    LAMA_LOG_DEBUG( logger, "benchmark is now shutdown" );
}

template<typename MatrixType>
CounterType PMetaSolverBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return mNumFloatingPointOperations;
}

template<typename MatrixType>
CounterType PMetaSolverBenchmark<MatrixType>::getProcessedBytes() const
{
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

#endif // LAMA_LAMAMETASOLVERBENCHMARK_HPP_
