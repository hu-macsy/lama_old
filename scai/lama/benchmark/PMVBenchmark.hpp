/**
 * @file PMVBenchmark.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Benchmark for parallel matrix-vector multiplication with different settings.
 * @author Thomas Brandes, Jiri Kraus
 * @date 11.05.2011
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
#include <scai/common/OpenMP.hpp>

#include <sstream>
#include <map>
#include <string>

// Benchmark: here we allow using namespace in .hpp file

using namespace scai;

using std::istringstream;
using std::map;
using std::string;

template<typename MatrixType>
class PMVBenchmark: 
 
    public  LAMAMPIBenchmark,
    private bf::Benchmark::Register< PMVBenchmark<MatrixType> >
{
public:

    typedef typename MatrixType::MatrixValueType ValueType;

    PMVBenchmark();

    PMVBenchmark( const std::string& arguments );

    PMVBenchmark( const PMVBenchmark<MatrixType>& other );

    virtual ~PMVBenchmark();

    virtual bf::Benchmark* copy() const;

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    /** Implementation of pure method Benchmark::getId()   */

    virtual const std::string& getId() const;

    static std::string createValue();

    static Benchmark* create()
    {
        return new PMVBenchmark();
    }

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

    static const std::string& getGroupId();

    common::unique_ptr<MatrixType> mMatrixA;
    common::unique_ptr<lama::DenseVector<ValueType> > mVectorX;
    common::unique_ptr<lama::DenseVector<ValueType> > mVectorY;

    dmemo::DistributionPtr mDistribution;
    hmemo::ContextPtr mContext;
    lama::Matrix::SyncKind mCommunicationKind;

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedBytesFloat;
    CounterType mNumProcessedBytesDouble;
};

template<typename MatrixType>
const std::string& PMVBenchmark<MatrixType>::getGroupId()
{
    static const std::string id = "PMV";
    return id;
}

template<typename MatrixType>
PMVBenchmark<MatrixType>::PMVBenchmark() :

    LAMAMPIBenchmark( getId(), getGroupId() ),
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
    static std::string id = createValue();
    return id;
}

template<typename MatrixType>
std::string PMVBenchmark<MatrixType>::createValue() 
{
    std::string value ( "PMV_" );
    value += MatrixType::typeName();
    return value;
}

template<typename MatrixType>
void PMVBenchmark<MatrixType>::initialize()
{
    //device initialization + comunicator

    map<string, string> tokens;

    getConfig( tokens );

    int noThreads = 1;
    int devNo = -1; // default device

    if ( tokens.count( "CUDA" ) > 0 )
    {
        istringstream tokenStream( tokens["CUDA"] );
        tokenStream >> devNo;
    }

    std::ostringstream idStream;
    idStream << " (Proc " << mComm->getRank() << " LOCAL=";

    if ( tokens["LOCAL"] == "CUDA" )
    {
        mContext = hmemo::Context::getContextPtr( hmemo::Context::CUDA, devNo );
        idStream << "CUDA:";
    }
    else
    {
        mContext = hmemo::Context::getContextPtr( hmemo::Context::Host );
        idStream << "HOST:";
    }

    if ( tokens.count( "THREADS" ) > 0 )
    {
        istringstream tokenStream( tokens["THREADS"] );
        tokenStream >> noThreads;
    }
    else
    {
        noThreads = 1;
    }

    idStream << "THREADS=" << noThreads << ":";

    if ( mContext->getType() == hmemo::Context::CUDA  )
    {
        idStream << "DEVICE=" << devNo << ":";
    }

    idStream << "COMM=";

    if ( tokens["COMM"] == "SYNC" )
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

    if ( tokens.count( "W" ) > 0 )
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

    const lama::LAMAInputSet& inputSet = bf::InputSetRegistry<lama::LAMAInputSet>::getRegistry().get( mInputSetId );

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
    const lama::LAMAInputSet& inputSet = bf::InputSetRegistry<lama::LAMAInputSet>::getRegistry().get( mInputSetId );

    //TODO: Complexity Calculation needs to be ported
    if ( ( typeid( MatrixType ) == typeid( lama::DenseMatrix<float> ) )
            || ( typeid( MatrixType ) == typeid( lama::DenseMatrix<double> ) ) )
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

    if ( typeid( ValueType ) == typeid( double ) )
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

        for ( int i = 0; i < correctResult.size(); i++ )
        {
            ValueType inputValue = lama::cast<ValueType>( ( *mVectorX )( i ) );
            ValueType correctValue = lama::cast<ValueType>( correctResult( i ) );
            ValueType computedValue = lama::cast<ValueType>( computedResult( i ) );
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
    catch ( bf::BFException& e )
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

    if ( sizeof( ValueType ) == sizeof( float ) )
    {
        return mNumProcessedBytesFloat;
    }
    else if ( sizeof( ValueType ) == sizeof( double ) )
    {
        return mNumProcessedBytesDouble;
    }

    return 0;
}
