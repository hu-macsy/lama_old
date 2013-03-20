/**
 * @file PCGBenchmark.hpp
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
 * @brief Benchmark for parallel CG implementation with different settings.
 * @author Jiri Kraus
 * @date 25.08.2011
 * $Id$
 */
#ifndef LAMA_LAMACGBENCHMARK_HPP_
#define LAMA_LAMACGBENCHMARK_HPP_

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

#include <lama/norm/MaxNorm.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/SpecializedJacobi.hpp>

#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>

#include <lama/solver/logger/OpenMPTimer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#ifdef LAMA_BUILD_CUDA
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <bench/LAMAMPIBenchmark.hpp>
#include <bench/LAMAInputSet.hpp>
#include <bench/LAMAInputSetComplexityVisitor.hpp>

using lama::LAMAInterfaceRegistry;
using lama::Context;
using lama::CG;
using lama::IterationCount;
using lama::ResidualThreshold;
using lama::Timer;
using lama::OpenMPTimer;
using lama::CommonLogger;
using lama::Logger;
using lama::LoggerPtr;
using lama::SpecializedJacobi;
using lama::CriterionPtr;
using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PCGBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PCGBenchmark();

    PCGBenchmark( const std::string& arguments );

    PCGBenchmark( const PCGBenchmark<MatrixType>& other );

    virtual ~PCGBenchmark();

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

    boost::shared_ptr<MatrixType> mMatrixA;

    typedef lama::DenseVector<ValueType> tVector;

    boost::shared_ptr<tVector> mSolution;
    boost::shared_ptr<tVector> mRhs;

    boost::shared_ptr<lama::CG> mCGSolver;

    lama::SolverPtr mPreconditioner;

    int mNumIterations;

    lama::DistributionPtr mDistribution;
    lama::ContextPtr mLocalContext;
    lama::ContextPtr mHaloContext;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
};

template<typename MatrixType>
const std::string& PCGBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PCGBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PCGBenchmark<MatrixType>::PCGBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mLocalContext(
        lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PCGBenchmark<MatrixType>::PCGBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mLocalContext(
        lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
            lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                0 ), mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PCGBenchmark<MatrixType>::PCGBenchmark( const PCGBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mLocalContext( lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
        lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
            0 ), mNumProcessedBytesDouble( 0 )
{

}

template<typename MatrixType>
PCGBenchmark<MatrixType>::~PCGBenchmark()
{
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PCGBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PCGBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PCGBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PCGBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PCGBenchmark<MatrixType>::getId() const
{
    return id();
}

template<typename MatrixType>
void PCGBenchmark<MatrixType>::initialize()
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

    mMatrixA.reset( new MatrixType( inputA ) );

    mMatrixA->redistribute( mDistribution, mDistribution );

    mMatrixA->setCommunicationKind( mCommunicationKind );
    mMatrixA->setContext( mLocalContext, mHaloContext );

    // redistribute vector Y of input set

    mRhs.reset( new lama::DenseVector<ValueType>( inputSet.getY(), mDistribution ) );

    mSolution.reset( new lama::DenseVector<ValueType>( mDistribution, 0.0 ) );

    mRhs->setContext( mLocalContext );
    mSolution->setContext( mLocalContext );

    const bool logging = false;

    if( logging )
    {
        LoggerPtr slogger(
            new CommonLogger( "<CG>: ", lama::LogLevel::completeInformation,
                              lama::LoggerWriteBehaviour::toConsoleOnly,
                              std::auto_ptr<Timer>( new OpenMPTimer() ) ) );

        mCGSolver.reset( new CG( "CG solver", slogger ) );
    }
    else
    {
        mCGSolver.reset( new CG( "CG solver" ) );
    }

    mNumIterations = 10;

    CriterionPtr minIterations( new IterationCount( mNumIterations ) );

    mCGSolver->setStoppingCriterion( minIterations );

    {
        SpecializedJacobi* jacobi = new SpecializedJacobi( "Preconditioner" );

        CriterionPtr minJacobiIterations( new IterationCount( 2 ) );

        jacobi->setStoppingCriterion( minJacobiIterations );
        jacobi->setOmega( 0.5 );

        mPreconditioner.reset( jacobi );
    }

    if( mPreconditioner )
    {
        mCGSolver->setPreconditioner( mPreconditioner );
    }

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << mLocalContext << ", halo at " << mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PCGBenchmark<MatrixType>::setUp()
{
    mCGSolver->initialize( *mMatrixA );

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
void PCGBenchmark<MatrixType>::execute()
{
    LAMA_LOG_INFO( logger, "execute benchmark" );

    mCGSolver->solve( *mSolution, *mRhs );
}

template<typename MatrixType>
void PCGBenchmark<MatrixType>::tearDown()
{
    LAMA_LOG_INFO( logger, "tear down" );

    lama::ContextPtr host = lama::ContextFactory::getContext( Context::Host );

    mSolution->prefetch( host );
    mSolution->wait();

    LAMA_LOG_INFO( logger, "tearDown done, Y at Host" );
}

template<typename MatrixType>
void PCGBenchmark<MatrixType>::shutdown()
{
    LAMA_LOG_INFO( logger, "shutdown benchmark" );

    const LAMAInputSet& inputSet = bf::InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    mCGSolver.reset();
    mMatrixA.reset();
    mSolution.reset();
    mRhs.reset();

    //TODO: Complexity Calculation needs to be ported

    LAMAInputSetComplexityVisitor::getCGComplexity( inputSet.getA(), mNumFloatingPointOperations,
            mNumProcessedBytesFloat, mNumProcessedBytesDouble );

    LAMA_LOG_DEBUG( logger, "benchmark is now shutdown" );
}

template<typename MatrixType>
CounterType PCGBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return mNumIterations * mNumFloatingPointOperations;
}

template<typename MatrixType>
CounterType PCGBenchmark<MatrixType>::getProcessedBytes() const
{
    if( sizeof(ValueType) == sizeof(float) )
    {
        return mNumIterations * mNumProcessedBytesFloat;
    }
    else if( sizeof(ValueType) == sizeof(double) )
    {
        return mNumIterations * mNumProcessedBytesDouble;
    }
    return 0;
}

#endif // LAMA_LAMACGBENCHMARK_HPP_
