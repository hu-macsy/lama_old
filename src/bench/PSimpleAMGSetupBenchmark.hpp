/**
 * @file PSimpleAMGSetupBenchmark.hpp
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
 * @brief Benchmark for parallel AMG Setup with different settings.
 * @author Malte Foerster
 * @date 16.11.2011
 * $Id$
 */
#ifndef LAMA_LAMASIMPLEAMGSETUPBENCHMARK_HPP_
#define LAMA_LAMASIMPLEAMGSETUPBENCHMARK_HPP_

#include <bench/LAMAMPIBenchmark.hpp>
#include <bench/LAMAInputSet.hpp>
#include <bench/LAMAInputSetComplexityVisitor.hpp>

#include <lama/DenseVector.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>
#include <lama/DefaultHostContextManager.hpp>

#include <lama/distribution/GenBlockDistribution.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <lama/norm/MaxNorm.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/SimpleAMG.hpp>

#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/tracing.hpp>
#include <lama/solver/logger/OpenMPTimer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#ifdef LAMA_BUILD_CUDA
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <map>
#include <string>

using lama::LAMAInterfaceRegistry;
using lama::Context;
using lama::CG;
using lama::SimpleAMG;
using lama::CriterionPtr;
using lama::IterationCount;
using lama::ResidualThreshold;
using lama::Criterion;
using lama::L2Norm;
using lama::Timer;
using lama::OpenMPTimer;
using lama::CommonLogger;
using lama::Logger;
using lama::LoggerPtr;
using lama::Scalar;
using lama::Context;
using lama::CSRSparseMatrix;
using lama::ELLSparseMatrix;
using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PSimpleAMGSetupBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PSimpleAMGSetupBenchmark();

    PSimpleAMGSetupBenchmark( const std::string& arguments );

    PSimpleAMGSetupBenchmark( const PSimpleAMGSetupBenchmark<MatrixType>& other );

    virtual ~PSimpleAMGSetupBenchmark();

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

    lama::Matrix* mMatrixA;

    lama::SimpleAMG* mAMGSolver;

    int mNumIterations;
    int mMaxLevels;
    bool mSolverLogging;
    int mHostOnlyLevel;
    std::string mArguments;
    mutable std::string mMyId;

    lama::DistributionPtr mDistribution;
    lama::ContextPtr mLocalContext;
    lama::ContextPtr mHaloContext;
    lama::ContextType mSmootherContextType;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
    std::auto_ptr<L2Norm> mL2Norm;
};

template<typename MatrixType>
const std::string& PSimpleAMGSetupBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PSimpleAMGSetupBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PSimpleAMGSetupBenchmark<MatrixType>::PSimpleAMGSetupBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mMatrixA( 0 ), mAMGSolver(
        0 ), mNumIterations( 1 ), mMaxLevels( 25 ), mSolverLogging( false ), mHostOnlyLevel(
            std::numeric_limits<int>::max() ), mLocalContext(
                lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mSmootherContextType(
                        Context::MaxContext ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat( 0 ), mNumProcessedBytesDouble(
                            0 ), mL2Norm( new L2Norm() )
{
}

template<typename MatrixType>
PSimpleAMGSetupBenchmark<MatrixType>::PSimpleAMGSetupBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mMatrixA(
        0 ), mAMGSolver( 0 ), mNumIterations( 1 ), mMaxLevels( 25 ), mSolverLogging( false ), mHostOnlyLevel(
            std::numeric_limits<int>::max() ), mArguments( arguments ), mLocalContext(
                lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mSmootherContextType(
                        Context::MaxContext ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat( 0 ), mNumProcessedBytesDouble(
                            0 ), mL2Norm( new L2Norm() )
{
    std::string::size_type firstSep = arguments.find( "_" );
    std::string::size_type secondSep = arguments.find( "_", firstSep + 1 );
    std::string::size_type thirdSep = arguments.find( "_", secondSep + 1 );
    std::string arg1 = arguments.substr( 0, firstSep );
    std::string arg2 = arguments.substr( firstSep + 1, 1 );
    std::string arg3 = arguments.substr( secondSep + 1, 1 );
    std::string arg4 = arguments.substr( thirdSep + 1 );
    {
        std::istringstream arg( arg1 );
        arg >> mMaxLevels;
    }
    if( arg2 != "" )
    {
        std::istringstream arg( arg2 );
        arg >> mSolverLogging;
    }
    if( arg3 != "" )
    {
        char c;
        std::istringstream arg( arg3 );
        arg >> c;
        if( c == 'C' )
        {
            mSmootherContextType = Context::CUDA;
        }
        else if( c == 'H' )
        {
            mSmootherContextType = Context::Host;
        }
    }
    if( arg4 != "" )
    {
        std::istringstream arg( arg4 );
        arg >> mHostOnlyLevel;
    }
    mName += "(" + arguments + ")";
}

template<typename MatrixType>
PSimpleAMGSetupBenchmark<MatrixType>::PSimpleAMGSetupBenchmark( const PSimpleAMGSetupBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mMatrixA( 0 ), mAMGSolver( 0 ), mNumIterations( other.mNumIterations ), mMaxLevels(
        other.mMaxLevels ), mSolverLogging( other.mSolverLogging ), mHostOnlyLevel(
            other.mHostOnlyLevel ), mArguments( other.mArguments ), mLocalContext(
                other.mLocalContext ), mHaloContext( other.mHaloContext ), mSmootherContextType(
                    other.mSmootherContextType ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                        0 ), mNumProcessedBytesDouble( 0 ), mL2Norm( new L2Norm() )
{

}

template<typename MatrixType>
PSimpleAMGSetupBenchmark<MatrixType>::~PSimpleAMGSetupBenchmark()
{
    if( mMatrixA != 0 )
    {
        delete mMatrixA;
    }
    mMatrixA = 0;
    if( mAMGSolver != 0 )
    {
        delete mAMGSolver;
    }
    mAMGSolver = 0;
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PSimpleAMGSetupBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PSimpleAMGSetupBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PSimpleAMGSetupBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PSimpleAMGSetupBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PSimpleAMGSetupBenchmark<MatrixType>::getId() const
{
    mMyId = id() + "(" + mArguments + ")";
    return mMyId;
}

template<typename MatrixType>
void PSimpleAMGSetupBenchmark<MatrixType>::initialize()
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

    LAMA_LOG_INFO( logger, "get matrix from input set" );

    const lama::CSRSparseMatrix<double>& inputA = inputSet.getA();

    LAMA_LOG_INFO( logger, "input matrix A = " << inputA );

    // mDistribution = lama::DistributionPtr ( new lama::GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );

    mDistribution = inputA.getDistributionPtr();

    // redistribute input matrix, is local copy in case of same distribution

    LAMA_LOG_INFO( logger, "create copy to desired matrix format" );

    if( mSmootherContextType == Context::MaxContext )
    {
        mSmootherContextType = mLocalContext->getType();
    }

    mMatrixA = new MatrixType( inputA );

    LAMA_LOG_INFO( logger, "redistribute matrix / should not do anything" );
    mMatrixA->redistribute( mDistribution, mDistribution );

    mMatrixA->setCommunicationKind( mCommunicationKind );
    mMatrixA->setContext( mLocalContext, mHaloContext );

    LAMA_LOG_INFO( logger, "create AMG Solver" );

    if( mSolverLogging )
    {
        bool ignoreRankInLog = ( mDistribution->getCommunicatorPtr()->getRank() != 0 );

        LoggerPtr amglogger(
            new CommonLogger( "<SimpleAMG>: ", lama::LogLevel::advancedInformation,
                              lama::LoggerWriteBehaviour::toConsoleOnly,
                              std::auto_ptr<Timer>( new OpenMPTimer() ), ignoreRankInLog ) );

        mAMGSolver = new SimpleAMG( "SimpleAMG solver", amglogger );
    }
    else
    {
        mAMGSolver = new SimpleAMG( "SimpleAMG solver" );
    }

    lama::ContextPtr smootherContext = lama::ContextFactory::getContext( mSmootherContextType, devNo );

    mAMGSolver->setSmootherContext( smootherContext );
    mAMGSolver->setHostOnlyLevel( mHostOnlyLevel );

    mAMGSolver->setMaxLevels( mMaxLevels );

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << *mLocalContext << ", halo at " << *mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PSimpleAMGSetupBenchmark<MatrixType>::setUp()
{
    LAMA_LOG_INFO( logger, "enter Benchmark::setUp" );

    mMatrixA->prefetch(); // prefetch local storage to local and halo storage to halo context
    mMatrixA->wait();

    lama::DefaultHostContextManager::setAsCurrent();
}

template<typename MatrixType>
void PSimpleAMGSetupBenchmark<MatrixType>::execute()
{
    LAMA_REGION( "PSimpleAMGSetupBenchmark::execute" );
    LAMA_LOG_INFO( logger, "execute benchmark" );

    mAMGSolver->initialize( *mMatrixA );

    LAMA_LOG_INFO( logger, "Initialize done for AMG" );
}

template<typename MatrixType>
void PSimpleAMGSetupBenchmark<MatrixType>::tearDown()
{
    LAMA_LOG_INFO( logger, "tear down" );
}

template<typename MatrixType>
void PSimpleAMGSetupBenchmark<MatrixType>::shutdown()
{
    LAMA_LOG_INFO( logger, "shutdown benchmark" );

    if( mAMGSolver != 0 )
    {
        delete mAMGSolver;
    }
    mAMGSolver = 0;

    if( mMatrixA != 0 )
    {
        delete mMatrixA;
    }
    mMatrixA = 0;

    LAMA_LOG_DEBUG( logger, "benchmark is now shutdown" );
}

template<typename MatrixType>
CounterType PSimpleAMGSetupBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return 0;
}

template<typename MatrixType>
CounterType PSimpleAMGSetupBenchmark<MatrixType>::getProcessedBytes() const
{
    return 0;
}

#endif // LAMA_LAMASIMPLEAMGSETUPBENCHMARK_HPP_
