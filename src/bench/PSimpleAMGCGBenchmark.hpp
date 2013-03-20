/**
 * @file PSimpleAMGCGBenchmark.hpp
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
 * @brief Benchmark for parallel AMG preconditioned CG implementation
 * with different settings.
 * @author Malte Foerster
 * @date 16.11.2011
 * $Id$
 */
#ifndef LAMA_LAMASIMPLEAMGCGBENCHMARK_HPP_
#define LAMA_LAMASIMPLEAMGCGBENCHMARK_HPP_

#include <sstream>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>

#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>
#include <lama/tracing.hpp>

#include <lama/distribution/GenBlockDistribution.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <lama/norm/MaxNorm.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/SimpleAMG.hpp>

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
using lama::SimpleAMG;
using lama::CriterionPtr;
using lama::IterationCount;
using lama::ResidualThreshold;
using lama::Criterion;
using lama::NormPtr;
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
class PSimpleAMGCGBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PSimpleAMGCGBenchmark();

    PSimpleAMGCGBenchmark( const std::string& arguments );

    PSimpleAMGCGBenchmark( const PSimpleAMGCGBenchmark<MatrixType>& other );

    virtual ~PSimpleAMGCGBenchmark();

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

    lama::DenseVector<ValueType>* mSolution;
    lama::DenseVector<ValueType>* mRhs;

    lama::CG* mCGSolver;
    boost::shared_ptr<lama::SimpleAMG> mAMGSolver;

    int mNumIterations;
    int mMaxLevels;
    bool mSolverLogging;
    int mHostOnlyLevel;
    int mReplicatedLevel;
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
    lama::NormPtr mL2Norm;
};

template<typename MatrixType>
const std::string& PSimpleAMGCGBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PSimpleAMGCGBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PSimpleAMGCGBenchmark<MatrixType>::PSimpleAMGCGBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mMatrixA( 0 ), mSolution(
        0 ), mRhs( 0 ), mCGSolver( 0 ), mNumIterations( 1 ), mMaxLevels( 25 ), mSolverLogging(
            false ), mHostOnlyLevel( std::numeric_limits<int>::max() ), mReplicatedLevel(
                std::numeric_limits<int>::max() ), mLocalContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                        lama::ContextFactory::getContext( Context::Host ) ), mSmootherContextType(
                            Context::MaxContext ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat( 0 ), mNumProcessedBytesDouble(
                                0 ), mL2Norm( new L2Norm() )
{
}

template<typename MatrixType>
PSimpleAMGCGBenchmark<MatrixType>::PSimpleAMGCGBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mMatrixA(
        0 ), mSolution( 0 ), mRhs( 0 ), mCGSolver( 0 ), mNumIterations( 1 ), mMaxLevels( 25 ), mSolverLogging(
            false ), mHostOnlyLevel( std::numeric_limits<int>::max() ), mReplicatedLevel(
                std::numeric_limits<int>::max() ), mArguments( arguments ), mLocalContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                        lama::ContextFactory::getContext( Context::Host ) ), mSmootherContextType(
                            Context::MaxContext ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat( 0 ), mNumProcessedBytesDouble(
                                0 ), mL2Norm( new L2Norm() )
{
    std::string::size_type firstSep = arguments.find( "_" );
    std::string::size_type secondSep = std::string::npos;
    if( firstSep != std::string::npos )
    {
        secondSep = arguments.find( "_", firstSep + 1 );
    }
    std::string::size_type thirdSep = std::string::npos;
    if( secondSep != std::string::npos )
    {
        thirdSep = arguments.find( "_", secondSep + 1 );
    }
    std::string::size_type fourthSep = std::string::npos;
    if( thirdSep != std::string::npos )
    {
        fourthSep = arguments.find( "_", thirdSep + 1 );
    }
    std::string arg1 = arguments.substr( 0, firstSep );
    std::string arg2 = "";
    if( firstSep != std::string::npos )
    {
        arg2 = arguments.substr( firstSep + 1, 1 );
    }
    std::string arg3 = "";
    if( secondSep != std::string::npos )
    {
        arg3 = arguments.substr( secondSep + 1, 1 );
    }
    std::string arg4 = "";
    if( thirdSep != std::string::npos )
    {
        arg4 = arguments.substr( thirdSep + 1, fourthSep - thirdSep );
    }
    std::string arg5 = "";
    if( fourthSep != std::string::npos )
    {
        arg5 = arguments.substr( fourthSep + 1, 100 );
    }
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
    if( arg5 != "" )
    {
        std::istringstream arg( arg5 );
        arg >> mReplicatedLevel;
    }

    LAMA_LOG_INFO( logger,
                   "args: mMaxLevels = " << mMaxLevels << ", mSolverLogging = " << mSolverLogging << ", mSmootherContextType = " << mSmootherContextType << ", mHostOnlyLevel = " << mHostOnlyLevel << ", mReplicatedLevel = " << mReplicatedLevel );

    mName += "(" + arguments + ")";
}

template<typename MatrixType>
PSimpleAMGCGBenchmark<MatrixType>::PSimpleAMGCGBenchmark( const PSimpleAMGCGBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mMatrixA( 0 ), mSolution( 0 ), mRhs( 0 ), mCGSolver( 0 ), mNumIterations(
        other.mNumIterations ), mMaxLevels( other.mMaxLevels ), mSolverLogging(
            other.mSolverLogging ), mHostOnlyLevel( other.mHostOnlyLevel ), mReplicatedLevel(
                other.mReplicatedLevel ), mArguments( other.mArguments ), mLocalContext(
                    other.mLocalContext ), mHaloContext( other.mHaloContext ), mSmootherContextType(
                        other.mSmootherContextType ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                            0 ), mNumProcessedBytesDouble( 0 ), mL2Norm( new L2Norm() )
{

}

template<typename MatrixType>
PSimpleAMGCGBenchmark<MatrixType>::~PSimpleAMGCGBenchmark()
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
    if( mCGSolver != 0 )
    {
        delete mCGSolver;
    }
    mCGSolver = 0;
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PSimpleAMGCGBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PSimpleAMGCGBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PSimpleAMGCGBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PSimpleAMGCGBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PSimpleAMGCGBenchmark<MatrixType>::getId() const
{
    mMyId = id() + "(" + mArguments + ")";
    return mMyId;
}

template<typename MatrixType>
void PSimpleAMGCGBenchmark<MatrixType>::initialize()
{
    LAMA_REGION( "initialize" );

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

    mDistribution = inputA.getDistributionPtr();

    int useWeightedDistribution = 0;
    if( weight != 1.0 )
    {
        useWeightedDistribution = 1;
    }

    useWeightedDistribution = mDistribution->getCommunicator().max( useWeightedDistribution );

    if( useWeightedDistribution ) {
        mDistribution = lama::DistributionPtr( new lama::GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );
    }

    // redistribute input matrix, is local copy in case of same distribution

    LAMA_LOG_INFO( logger, "create copy to desired matrix format" );

    lama::ContextPtr smootherContext;

    if( mSmootherContextType != Context::MaxContext )
    {
        smootherContext = lama::ContextFactory::getContext( mSmootherContextType, devNo );
    }

    if( group() == LAMAInputSetComplexityVisitor::ELLCSRSIMPLEAMGCG )
    {
        if( mLocalContext->getType() == Context::CUDA )
        {
            mMatrixA = new lama::ELLSparseMatrix<ValueType>( inputA );
        }
        else
        {
            mMatrixA = new lama::CSRSparseMatrix<ValueType>( inputA );
        }
    }
    else
    {
        mMatrixA = new MatrixType( inputA );
    }

    if( !useWeightedDistribution )
    {
        LAMA_LOG_INFO( logger, "redistribute matrix / should not do anything" );
    }
    mMatrixA->redistribute( mDistribution, mDistribution );

    mMatrixA->setCommunicationKind( mCommunicationKind );
    mMatrixA->setContext( mLocalContext, mHaloContext );

    // redistribute vector Y of input set

    LAMA_LOG_INFO( logger, "create vectors for rhs / solution" );

    mRhs = new lama::DenseVector<ValueType>( inputSet.getY(), mDistribution );

    mSolution = new lama::DenseVector<ValueType>( mDistribution, 0.0 );

    mRhs->setContext( mLocalContext );
    mSolution->setContext( mLocalContext );

    LAMA_LOG_INFO( logger, "create CG & AMG Solver" );

    if( mSolverLogging )
    {
        bool ignoreRankInLog = ( mDistribution->getCommunicatorPtr()->getRank() != 0 );

        LoggerPtr cglogger(
            new CommonLogger( "<CG/SimpleAMG>: ", lama::LogLevel::convergenceHistory,
                              lama::LoggerWriteBehaviour::toConsoleOnly,
                              std::auto_ptr<Timer>( new OpenMPTimer() ), ignoreRankInLog ) );
        LoggerPtr amglogger(
            new CommonLogger( "<SimpleAMG>: ", lama::LogLevel::advancedInformation,
                              lama::LoggerWriteBehaviour::toConsoleOnly,
                              std::auto_ptr<Timer>( new OpenMPTimer() ), ignoreRankInLog ) );

        mCGSolver = new CG( "CG solver", cglogger );
        mAMGSolver.reset( new SimpleAMG( "SimpleAMG solver", amglogger ) );
    }
    else
    {
        mCGSolver = new CG( "CG solver" );
        mAMGSolver.reset( new SimpleAMG( "SimpleAMG solver" ) );
    }

    mAMGSolver->setSmootherContext( smootherContext );
    mAMGSolver->setHostOnlyLevel( mHostOnlyLevel );
    mAMGSolver->setReplicatedLevel( mReplicatedLevel );

    mNumIterations = 10;
    CriterionPtr minIterations( new IterationCount( mNumIterations ) );
    CriterionPtr resThreashold( new ResidualThreshold( mL2Norm, Scalar( 0.0 ), ResidualThreshold::Absolute ) );

    CriterionPtr cgCriterion = minIterations || resThreashold;
    mCGSolver->setStoppingCriterion( cgCriterion );

    mAMGSolver->setMaxLevels( mMaxLevels );

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << *mLocalContext << ", halo at " << *mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PSimpleAMGCGBenchmark<MatrixType>::setUp()
{
    LAMA_REGION( "setUp_PSimpleAMGCG" );

    LAMA_LOG_INFO( logger, "enter Benchmark::setUp" );

    mCGSolver->setPreconditioner( mAMGSolver );

    mCGSolver->initialize( *mMatrixA );

    LAMA_LOG_INFO( logger, "Initialize done for CG and its preconditioner (AMG)" );

    mAMGSolver->setLogLevel( lama::LogLevel::noLogging );

    mSolution->prefetch( mLocalContext );
    mRhs->prefetch( mLocalContext );
    mMatrixA->prefetch(); // prefetch local storage to local and halo storage to halo context
    mSolution->wait();
    mRhs->wait();
    mMatrixA->wait();

    LAMA_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : Solution, Rhs and  A local at " << *mLocalContext << " and A halo at " << *mHaloContext );
}

template<typename MatrixType>
void PSimpleAMGCGBenchmark<MatrixType>::execute()
{
    LAMA_REGION( "execute_PSimpleAMGCG" );

    LAMA_LOG_INFO( logger, "execute benchmark" );

    mCGSolver->solve( *mSolution, *mRhs );
}

template<typename MatrixType>
void PSimpleAMGCGBenchmark<MatrixType>::tearDown()
{
    LAMA_REGION( "tearDown_PSimpleAMGCG" );

    LAMA_LOG_INFO( logger, "tear down" );

    lama::ContextPtr host = lama::ContextFactory::getContext( Context::Host );

    mRhs->prefetch( host );
    mRhs->wait();

    LAMA_LOG_INFO( logger, "tearDown done, Y at Host" );
}

template<typename MatrixType>
void PSimpleAMGCGBenchmark<MatrixType>::shutdown()
{
    LAMA_REGION( "shutdown_PSimpleAMGCG" );

    LAMA_LOG_INFO( logger, "shutdown benchmark" );

    //Complexity of AMG solver
    int numLevels = mAMGSolver->getNumLevels();

    CounterType numFloatingPointOperations = 0;
    CounterType processedBytesFloat = 0;
    CounterType processedBytesDouble = 0;

    // Complexity of CG
    LAMAInputSetComplexityVisitor::getCGComplexity( mAMGSolver->getGalerkin( 0 ), numFloatingPointOperations,
            processedBytesFloat, processedBytesDouble );

    for( int i = 0; i < ( numLevels - 1 ); ++i )
    {
        CounterType numFlops = 0;
        CounterType bytesFloat = 0;
        CounterType bytesDouble = 0;

        const lama::Matrix& galerkin = mAMGSolver->getGalerkin( i );
        const lama::Matrix& restriction = mAMGSolver->getRestriction( i );
        const lama::Matrix& interpolation = mAMGSolver->getInterpolation( i );
        int nnuGalerkin = galerkin.getNumRows();
        int nnuRestriction = restriction.getNumRows();
        int nnuInterpolation = interpolation.getNumRows();

        //Complexity of Smoothing (2 pre and 2 post smoothing steps (2times2))
        //Jacobi has as similar complexity
        LAMAInputSetComplexityVisitor::getMVComplexity( galerkin, numFlops, bytesFloat, bytesDouble );
        //Add extra accesses to old solution and rhs
        bytesFloat += 2 * nnuGalerkin * sizeof(float);
        bytesDouble += 2 * nnuGalerkin * sizeof(double);
        //Flops of relaxation
        numFlops += nnuGalerkin * 4;
        //4 Times (2 pre and 2 post smoothing steps (2times2))
        numFloatingPointOperations += 4 * numFlops;
        processedBytesFloat += 4 * bytesFloat;
        processedBytesDouble += 4 * bytesDouble;

        //Complexity of residue calculation
        LAMAInputSetComplexityVisitor::getMVComplexity( galerkin, numFlops, bytesFloat, bytesDouble );
        //Calculation of the residue is v = y - A*x not SpMV (v = A*x)
        numFloatingPointOperations += ( numFlops + nnuGalerkin );
        processedBytesFloat += ( bytesFloat + nnuGalerkin * sizeof(float) );
        processedBytesDouble += ( bytesDouble + nnuGalerkin * sizeof(double) );
        //Complexity of Restriction
        LAMAInputSetComplexityVisitor::getMVComplexity( restriction, numFlops, bytesFloat, bytesDouble );
        numFloatingPointOperations += numFlops;
        processedBytesFloat += bytesFloat;
        processedBytesDouble += bytesDouble;

        //Complexity of setting coarse solution to 0.0
        processedBytesFloat += nnuRestriction * sizeof(float);
        processedBytesDouble += nnuRestriction * sizeof(double);

        //Complexity of interpolation
        LAMAInputSetComplexityVisitor::getMVComplexity( interpolation, numFlops, bytesFloat, bytesDouble );
        //Calculation of the interpolation is v = y + A*x not SpMV (v = A*x)
        numFloatingPointOperations += ( numFlops + nnuInterpolation );
        processedBytesFloat += ( bytesFloat + nnuInterpolation * sizeof(float) );
        processedBytesDouble += ( bytesDouble + nnuInterpolation * sizeof(double) );
    }
    //Complexity of coarse level solver
    {
        //As an approximation consider the inverse of Galerkin to be
        //full. The complexity of the coarse level solver is small anyway.
        const int nnu = mAMGSolver->getGalerkin( numLevels - 1 ).getNumRows();
        numFloatingPointOperations += 2 * nnu * nnu - nnu;
        processedBytesFloat += ( nnu * nnu + 2 * nnu ) * sizeof(float);
        processedBytesDouble += ( nnu * nnu + 2 * nnu ) * sizeof(double);
    }

    mNumFloatingPointOperations = numFloatingPointOperations;
    mNumProcessedBytesFloat = processedBytesFloat;
    mNumProcessedBytesDouble = processedBytesDouble;

    LAMA_LOG_DEBUG( logger, "benchmark is now shutdown" );

    if( mCGSolver != 0 )
    {
        delete mCGSolver;
    }
    mCGSolver = 0;

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
CounterType PSimpleAMGCGBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return mNumIterations * mNumFloatingPointOperations;
}

template<typename MatrixType>
CounterType PSimpleAMGCGBenchmark<MatrixType>::getProcessedBytes() const
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

#endif // LAMA_LAMASIMPLEAMGCGBENCHMARK_HPP_
