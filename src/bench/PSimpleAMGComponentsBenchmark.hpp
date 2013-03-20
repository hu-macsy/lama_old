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
 * @brief Benchmark for parallel AMG Components
 * with different settings.
 * @author Malte Foerster
 * @date 16.11.2011
 * $Id$
 */
#ifndef LAMA_LAMASIMPLEAMGCOMPONENTSBENCHMARK_HPP_
#define LAMA_LAMASIMPLEAMGCOMPONENTSBENCHMARK_HPP_

#include <bench/LAMAMPIBenchmark.hpp>
#include <bench/LAMAInputSet.hpp>
#include <bench/LAMAInputSetComplexityVisitor.hpp>

#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>

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
#include <lama/norm/L2Norm.hpp>

#include <lama/solver/SimpleAMG.hpp>

#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>

#include <lama/solver/logger/OpenMPTimer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#include <lama/tracing.hpp>

#ifdef LAMA_BUILD_CUDA
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <map>
#include <string>

using lama::LAMAInterfaceRegistry;
using lama::Context;
using lama::SimpleAMG;
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
class PSimpleAMGComponentsBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    PSimpleAMGComponentsBenchmark();

    PSimpleAMGComponentsBenchmark( const std::string& arguments );

    PSimpleAMGComponentsBenchmark( const PSimpleAMGComponentsBenchmark<MatrixType>& other );

    virtual ~PSimpleAMGComponentsBenchmark();

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

    lama::SimpleAMG* mAMGSolver;

    int mMaxLevel;
    int mCheckLevel;
    bool mSolverLogging;

    bool mRunInterpolation;
    bool mRunRestriction;
    bool mRunJacobi;

    std::string mArguments;
    mutable std::string mMyId;

    lama::DistributionPtr mDistribution;
    lama::ContextPtr mLocalContext;
    lama::ContextPtr mHaloContext;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
    std::auto_ptr<L2Norm> mL2Norm;
};

template<typename MatrixType>
const std::string& PSimpleAMGComponentsBenchmark<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& PSimpleAMGComponentsBenchmark<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
PSimpleAMGComponentsBenchmark<MatrixType>::PSimpleAMGComponentsBenchmark()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), mMatrixA( 0 ), mSolution(
        0 ), mRhs( 0 ), mAMGSolver( 0 ), mMaxLevel( 25 ), mCheckLevel( 25 ), mSolverLogging(
            false ), mRunInterpolation( false ), mRunRestriction( false ), mRunJacobi( false ), mLocalContext(
                lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                        0 ), mNumProcessedBytesDouble( 0 ), mL2Norm( new L2Norm() )
{
}

template<typename MatrixType>
PSimpleAMGComponentsBenchmark<MatrixType>::PSimpleAMGComponentsBenchmark( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ), mMatrixA(
        0 ), mSolution( 0 ), mRhs( 0 ), mAMGSolver( 0 ), mMaxLevel( 25 ), mCheckLevel( 25 ), mSolverLogging(
            false ), mRunInterpolation( false ), mRunRestriction( false ), mRunJacobi( false ), mArguments(
                arguments ), mLocalContext( lama::ContextFactory::getContext( Context::Host ) ), mHaloContext(
                    lama::ContextFactory::getContext( Context::Host ) ), mNumFloatingPointOperations( 0 ), mNumProcessedBytesFloat(
                        0 ), mNumProcessedBytesDouble( 0 ), mL2Norm( new L2Norm() )
{
    std::string::size_type firstSep = arguments.find( "_" );
    std::string::size_type secondSep = arguments.find( "_", firstSep + 1 );
    std::string::size_type thirdSep = arguments.find( "_", secondSep + 1 );
    std::string arg1 = arguments.substr( 0, firstSep );
    std::string arg2 = arguments.substr( firstSep + 1, secondSep );
    std::string arg3 = arguments.substr( secondSep + 1, 1 );
    std::string arg4 = arguments.substr( thirdSep + 1, 1 );
    {
        std::istringstream arg( arg1 );
        arg >> mMaxLevel;
    }
    {
        std::istringstream arg( arg2 );
        arg >> mCheckLevel;
    }
    {
        std::istringstream arg( arg3 );
        arg >> mSolverLogging;
    }
    {
        char c;
        std::istringstream arg( arg4 );
        arg >> c;
        if( c == 'X' )
        {
            mRunInterpolation = true;
            mRunRestriction = true;
            mRunJacobi = true;
        }
        else if( c == 'I' )
        {
            mRunInterpolation = true;
        }
        else if( c == 'R' )
        {
            mRunRestriction = true;
        }
        else if( c == 'J' )
        {
            mRunJacobi = true;
        }
    }

    // in case check coarse grid solver
    if( mMaxLevel == mCheckLevel + 1 )
    {
        mRunInterpolation = false;
        mRunRestriction = false;
        mRunJacobi = false;
    }

    LAMA_LOG_INFO( logger,
                   "args: mMaxLevel = " << mMaxLevel << ", mCheckLevel = " << mCheckLevel << ", mSolverLogging = " << mSolverLogging << ", mRunInterpolation = " << mRunInterpolation << ", mRunRestriction = " << mRunRestriction << ", mRunJacobi = " << mRunJacobi );

    mName += "(" + arguments + ")";
}

template<typename MatrixType>
PSimpleAMGComponentsBenchmark<MatrixType>::PSimpleAMGComponentsBenchmark(
    const PSimpleAMGComponentsBenchmark<MatrixType>& other )
    : LAMAMPIBenchmark( other ), mMatrixA( 0 ), mSolution( 0 ), mRhs( 0 ), mAMGSolver( 0 ), mMaxLevel(
        other.mMaxLevel ), mCheckLevel( other.mCheckLevel ), mSolverLogging(
            other.mSolverLogging ), mRunInterpolation( other.mRunInterpolation ), mRunRestriction(
                other.mRunRestriction ), mRunJacobi( other.mRunJacobi ), mArguments( other.mArguments ), mLocalContext(
                    other.mLocalContext ), mHaloContext( other.mHaloContext ), mNumFloatingPointOperations(
                        0 ), mNumProcessedBytesFloat( 0 ), mNumProcessedBytesDouble( 0 ), mL2Norm(
                            new L2Norm() )
{

}

template<typename MatrixType>
PSimpleAMGComponentsBenchmark<MatrixType>::~PSimpleAMGComponentsBenchmark()
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
    if( mAMGSolver != 0 )
    {
        delete mAMGSolver;
    }
    mAMGSolver = 0;
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> PSimpleAMGComponentsBenchmark<MatrixType>::copy() const
{
    bf::Benchmark* b = new PSimpleAMGComponentsBenchmark<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short PSimpleAMGComponentsBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool PSimpleAMGComponentsBenchmark<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& PSimpleAMGComponentsBenchmark<MatrixType>::getId() const
{
    mMyId = id() + "(" + mArguments + ")";
    return mMyId;
}

template<typename MatrixType>
void PSimpleAMGComponentsBenchmark<MatrixType>::initialize()
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

    mMatrixA = new MatrixType( inputA );

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

    mAMGSolver->setSmootherContext( mLocalContext );

    // we need at least one more level than mCheckLevel no have a smoother
    // additionally increase by one cause we start counting at level 0
    if( mMaxLevel < mCheckLevel + 1 )
    {
        LAMA_LOG_ERROR( logger, "Level Hierarchy too small" );
    }
    mAMGSolver->setMaxLevels( mMaxLevel );

    LAMA_LOG_INFO( logger,
                   "Matrix: local at " << *mLocalContext << ", halo at " << *mHaloContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PSimpleAMGComponentsBenchmark<MatrixType>::setUp()
{
    LAMA_REGION( "setUp_PSimpleAMGComponents" );

    LAMA_LOG_INFO( logger, "enter Benchmark::setUp" );

    mAMGSolver->initialize( *mMatrixA );

    LAMA_LOG_INFO( logger, "Initialize done for AMG" );

    mAMGSolver->setLogLevel( lama::LogLevel::noLogging );

    mSolution->prefetch( mLocalContext );
    mRhs->prefetch( mLocalContext );
    mMatrixA->prefetch(); // prefetch local storage to local and halo storage to halo context
    mSolution->wait();
    mRhs->wait();
    mMatrixA->wait();

    // initialize vectors
    if( mCheckLevel != 0 )
    {
        lama::Vector& curSolution = mAMGSolver->getSolutionVector( mCheckLevel );
        lama::Vector& curRhs = mAMGSolver->getRhsVector( mCheckLevel );
        curSolution = 1.0;
        curRhs = 1.0;
    }

    // don't do this in case of InverseSolvertest
    if( mCheckLevel < mMaxLevel - 1 )
    {
        lama::Vector& curCoarseSolution = mAMGSolver->getSolutionVector( mCheckLevel + 1 );
        lama::Vector& curCoarseRhs = mAMGSolver->getRhsVector( mCheckLevel + 1 );
        curCoarseRhs = 1.0;
        curCoarseSolution = 1.0;
    }

    LAMA_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : Solution, Rhs and  A local at " << *mLocalContext << " and A halo at " << *mHaloContext );
}

template<typename MatrixType>
void PSimpleAMGComponentsBenchmark<MatrixType>::execute()
{
    LAMA_REGION( "execute_PSimpleAMGComponents" );

    const lama::Vector* curRhsPtr = mRhs;
    lama::Vector* curSolutionPtr = mSolution;
    if( mCheckLevel != 0 )
    {
        curSolutionPtr = &( mAMGSolver->getSolutionVector( mCheckLevel ) );
        curRhsPtr = &( mAMGSolver->getRhsVector( mCheckLevel ) );
    }
    if( mCheckLevel == mMaxLevel - 1 )
    {
        lama::Vector& curSolution = ( *curSolutionPtr );
        const lama::Vector& curRhs = ( *curRhsPtr );
        mAMGSolver->getCoarseLevelSolver().solve( curSolution, curRhs );
    }
    else
    {
        lama::Vector& curSolution = ( *curSolutionPtr );
        const lama::Vector& curRhs = ( *curRhsPtr );
        lama::Vector& curCoarseSolution = mAMGSolver->getSolutionVector( mCheckLevel + 1 );
        lama::Vector& curCoarseRhs = mAMGSolver->getRhsVector( mCheckLevel + 1 );

        LAMA_LOG_INFO( logger, "execute benchmark" );
        if( mRunInterpolation )
        {
            LAMA_LOG_INFO( logger, "execute interpolation" );
            const lama::Matrix& interpolation = mAMGSolver->getInterpolation( mCheckLevel );
            curSolution = curSolution + interpolation * curCoarseSolution;
        }
        if( mRunRestriction )
        {
            LAMA_LOG_INFO( logger, "execute restriction" );
            const lama::Matrix& restriction = mAMGSolver->getRestriction( mCheckLevel );
            curCoarseRhs = restriction * curSolution;
        }
        if( mRunJacobi )
        {
            LAMA_LOG_INFO( logger, "execute jacobi" );
            mAMGSolver->getSmoother( mCheckLevel ).solve( curSolution, curRhs );
        }
    }
}

template<typename MatrixType>
void PSimpleAMGComponentsBenchmark<MatrixType>::tearDown()
{
    LAMA_REGION( "tearDown_PSimpleAMGComponents" );

    LAMA_LOG_INFO( logger, "tear down" );

    lama::ContextPtr host = lama::ContextFactory::getContext( Context::Host );

    mRhs->prefetch( host );
    mRhs->wait();

    LAMA_LOG_INFO( logger, "tearDown done, Y at Host" );
}

template<typename MatrixType>
void PSimpleAMGComponentsBenchmark<MatrixType>::shutdown()
{
    LAMA_REGION( "shutdown_PSimpleAMGComponents" );

    LAMA_LOG_INFO( logger, "shutdown benchmark" );

    CounterType numFloatingPointOperations = 0;
    CounterType processedBytesFloat = 0;
    CounterType processedBytesDouble = 0;

    int i = mCheckLevel;

    CounterType numFlops = 0;
    CounterType bytesFloat = 0;
    CounterType bytesDouble = 0;

    if( mRunJacobi )
    {
        //Complexity of Smoothing (2 smoothing steps)
        const lama::Matrix& galerkin = mAMGSolver->getGalerkin( i );
        LAMAInputSetComplexityVisitor::getDefaultJacobiComplexity( 2, galerkin, numFlops, bytesFloat, bytesDouble );
        numFloatingPointOperations += numFlops;
        processedBytesFloat += bytesFloat;
        processedBytesDouble += bytesDouble;
    }

    if( mRunRestriction )
    {
        //Complexity of Restriction
        const lama::Matrix& restriction = mAMGSolver->getRestriction( i );
        LAMAInputSetComplexityVisitor::getMVComplexity( restriction, numFlops, bytesFloat, bytesDouble );
        numFloatingPointOperations += numFlops;
        processedBytesFloat += bytesFloat;
        processedBytesDouble += bytesDouble;
    }

    if( mRunInterpolation )
    {
        //Complexity of interpolation
        const lama::Matrix& interpolation = mAMGSolver->getInterpolation( i );
        int nnuInterpolation = interpolation.getNumRows();
        LAMAInputSetComplexityVisitor::getMVComplexity( interpolation, numFlops, bytesFloat, bytesDouble );
        //Calculation of the interpolation is v = y + A*x not SpMV (v = A*x)
        numFloatingPointOperations += ( numFlops + nnuInterpolation );
        processedBytesFloat += ( bytesFloat + nnuInterpolation * sizeof(float) );
        processedBytesDouble += ( bytesDouble + nnuInterpolation * sizeof(double) );
    }

    //Complexity of coarse level solver
    if( mCheckLevel == mMaxLevel - 1 )
    {
        //As an approximation consider the inverse of Galerkin to be
        //full. The complexity of the coarse level solver is small anyway.
        const int nnu = mAMGSolver->getGalerkin( mCheckLevel ).getNumRows();
        numFloatingPointOperations += 2 * nnu * nnu - nnu;
        processedBytesFloat += ( nnu * nnu + 2 * nnu ) * sizeof(float);
        processedBytesDouble += ( nnu * nnu + 2 * nnu ) * sizeof(double);
    }

    mNumFloatingPointOperations = numFloatingPointOperations;
    mNumProcessedBytesFloat = processedBytesFloat;
    mNumProcessedBytesDouble = processedBytesDouble;

    LAMA_LOG_DEBUG( logger, "benchmark is now shutdown" );

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
CounterType PSimpleAMGComponentsBenchmark<MatrixType>::getNumFloatingPointOperations() const
{
    return mNumFloatingPointOperations;
}

template<typename MatrixType>
CounterType PSimpleAMGComponentsBenchmark<MatrixType>::getProcessedBytes() const
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

#endif // LAMA_LAMASIMPLEAMGCOMPONENTSBENCHMARK_HPP_
