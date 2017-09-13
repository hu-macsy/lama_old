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

#pragma once

#include <scai/lama/benchmark/LAMAMPIBenchmark.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>
#include <scai/lama/benchmark/LAMAInputSetComplexityVisitor.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/dmemo/GenBlockDistribution.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <map>
#include <string>

// Benchmark: here we allow using namespace in .hpp file 

using namespace scai;

using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PMVBenchmark: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::MatrixValueType ValueType;

    static const std::string& id();

    PMVBenchmark();

    PMVBenchmark( const std::string& arguments );

    PMVBenchmark( const PMVBenchmark<MatrixType>& other );

    virtual ~PMVBenchmark();

    virtual bf::Benchmark* copy() const;

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

    common::unique_ptr<MatrixType> mMatrixA;
    common::unique_ptr<lama::DenseVector<ValueType> > mVectorX;
    common::unique_ptr<lama::DenseVector<ValueType> > mVectorY;

    dmemo::DistributionPtr mDistribution;
    hmemo::ContextPtr mContext;
    scai::lama::Matrix::SyncKind mCommunicationKind;

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
PMVBenchmark<MatrixType>::PMVBenchmark() : 

    LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) ), 
    mContext( hmemo::Context::getContextPtr( hmemo::Context::Host ) ),
    mNumFloatingPointOperations( 0 ), 
    mNumProcessedBytesFloat( 0 ), 
    mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark( const std::string& arguments ) : 

    LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments ),
    mContext( hmemo::Context::getContextPtr( hmemo::Context::Host ) ),
    mNumFloatingPointOperations( 0 ), 
    mNumProcessedBytesFloat( 0 ), 
    mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark( const PMVBenchmark<MatrixType>& other ) :

    LAMAMPIBenchmark( other ),
    mContext( hmemo::Context::getContextPtr( hmemo::Context::Host ) ),
    mNumFloatingPointOperations( 0 ), 
    mNumProcessedBytesFloat( 0 ), 
    mNumProcessedBytesDouble( 0 )
{
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::~PMVBenchmark()
{
    // Note: mMatrixA, mVectorX, mVectorY will be freed
}

template<typename MatrixType>
bf::Benchmark* PMVBenchmark<MatrixType>::copy() const
{
    return new PMVBenchmark<MatrixType>( *this );
}

template<typename MatrixType>
short PMVBenchmark<MatrixType>::getValueTypeSize() const
{
    return sizeof( ValueType );
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
        mContext = hmemo::Context::getContextPtr( hmemo::Context::CUDA, devNo );
        idStream << "CUDA:";
    }
    else
    {
        mContext = hmemo::Context::getContextPtr( hmemo::Context::Host );
        idStream << "HOST:";
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
    if( mContext->getType() == hmemo::Context::CUDA  )
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

    SCAI_LOG_INFO( logger, "get input set " << mInputSetId );

    const LAMAInputSet& inputSet = bf::InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    const lama::DenseVector<double>& inputX = inputSet.getX();
    const lama::CSRSparseMatrix<double>& inputA = inputSet.getA();

    SCAI_LOG_INFO( logger, "input matrix A = " << inputA );

    // mDistribution = lama::DistributionPtr ( new lama::GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );

    mDistribution = inputA.getRowDistributionPtr();

    SCAI_LOG_INFO( logger,
                   "General block distribution of matrix with weight " << weight << ", local size = " << mDistribution->getLocalSize() );

    mVectorX.reset( new lama::DenseVector<ValueType>( inputX, mDistribution ) );
    mVectorY.reset( new lama::DenseVector<ValueType>( mDistribution, 0.0 ) );

    SCAI_LOG_DEBUG( logger, "convert and redistribute the input matrix " << inputA );

    mMatrixA.reset( new MatrixType( inputA ) );

    mMatrixA->redistribute( mDistribution, mDistribution );

    SCAI_LOG_DEBUG( logger, "matrix " << *mMatrixA << " now available for MV" );

    mMatrixA->setContextPtr( mContext );
    mMatrixA->setCommunicationKind( mCommunicationKind );

    SCAI_LOG_INFO( logger,
                   "Matrix: context at " << mContext << ", comm = " << mCommunicationKind );
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::setUp()
{
    mVectorX->prefetch( mContext );
    mMatrixA->prefetch();
    mVectorX->wait();
    mMatrixA->wait();

    SCAI_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : X, A at " << mContext );
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::execute()
{
    lama::Vector& result = *mVectorY;
    const lama::Matrix& A = *mMatrixA;
    const lama::Vector& x = *mVectorX;
    result = A * x;
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::tearDown()
{
    hmemo::ContextPtr host = hmemo::Context::getHostPtr();
    mVectorY->prefetch( host );
    mVectorY->wait();
    SCAI_LOG_INFO( logger, "tearDown done, Y at Host" );
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

        SCAI_LOG_INFO( logger, "computed result: " << computedResult );

        // Note: there might be type conversion from double to float

        lama::DenseVector<ValueType> correctResult( result );

        SCAI_LOG_INFO( logger, "result ( by inputSet ): " << result );
        SCAI_LOG_INFO( logger, "correct result: " << correctResult );

#if defined( SCAI_LOG_TRACE_ENABLED)
        for ( int i = 0; i < correctResult.size(); i++)
        {
            ValueType inputValue = lama::cast<ValueType>( ( *mVectorX )( i ) );
            ValueType correctValue = lama::cast<ValueType>( correctResult(i) );
            ValueType computedValue = lama::cast<ValueType>( computedResult(i) );
            SCAI_LOG_TRACE( logger, i << ": correct = " << correctValue
                            << ", computed = " << computedValue
                            << ", input = " << inputValue );
        }
#endif

        lama::Vector& diff = correctResult;

        diff = computedResult - correctResult;

        lama::Scalar diffNorm = lama::maxNorm( diff );

        SCAI_LOG_INFO( logger, "max diff = " << diffNorm );

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
    typedef typename MatrixType::MatrixValueType ValueType;

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
