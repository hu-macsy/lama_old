/**
 * @file PMVBenchmark.cpp
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
 * @brief Implementation of methods for parallel benchmark MatrixTimesVector
 * @author Thomas Brandes
 * @date 21.09.2017
 */

#include <scai/lama/benchmark/PMVBenchmark.hpp>

#include <scai/common/macros/loop.hpp>

namespace scai
{

using std::istringstream;
using std::map;
using std::string;

namespace lama
{

const std::string& PMVBenchmark::getGroupId()
{
    static const std::string id = "PMV";
    return id;
}

PMVBenchmark::PMVBenchmark( const std::string& argument ) :

    LAMAMPIBenchmark( getCreateId(), getGroupId() ),
    mContext( hmemo::Context::getContextPtr( hmemo::Context::Host ) ),
    mNumFloatingPointOperations( 0 ),
    mNumProcessedBytesFloat( 0 ),
    mNumProcessedBytesDouble( 0 )
{
    // ToDo: args = CSR , double -> use as key to create matrix storage

    std::vector<std::string> argTokens;

    if ( argument == "" )
    {
        tokenize( argTokens, "CSR double", " ,:" );
    }
    else
    {
        tokenize( argTokens, argument, " ,:" );
    }

    SCAI_ASSERT_EQ_ERROR( argTokens.size(), 2, "format, type expected as argument" )

    Matrix::MatrixStorageFormat format = str2Format( argTokens[0].c_str() );

    common::scalar::ScalarType type = common::str2ScalarType( argTokens[1].c_str() );

    // allocate source and target storage of the required type

    mMatrixA.reset( Matrix::create( MatrixCreateKeyType( format, type ) ) );
    mVectorX.reset( _DenseVector::create( type ) );
    mVectorY.reset( _DenseVector::create( type ) );

    mArgument = argTokens[0] + ", " + argTokens[1];
}

PMVBenchmark::~PMVBenchmark()
{
    // Note: mMatrixA, mVectorX, mVectorY will be freed
}

short PMVBenchmark::getValueTypeSize() const
{
    return common::typeSize( mMatrixA->getValueType() );
}

bool PMVBenchmark::isThreadded() const
{
    return true;
}

const std::string& PMVBenchmark::getCreateId() const
{
    static std::string id = createValue();
    return id;
}

const std::string& PMVBenchmark::getArgument() const
{
    return mArgument;
}

std::string PMVBenchmark::createValue() 
{
    std::string value ( "PMatTimesVector" );
    return value;
}

void PMVBenchmark::initialize()
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
        mCommunicationKind = Matrix::SYNCHRONOUS;
        idStream << "SYNC:";
    }
    else
    {
        mCommunicationKind = Matrix::ASYNCHRONOUS;
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

    SCAI_LOG_INFO( logger, "get input set by this id " << mInputSetId );

    // mInputSetId = "inputSetId( argument )" must also be parsed ( "inputSetId", "argument" )

    mInputSet.reset( benchmark::InputSet::createWithArgument( mInputSetId ) );

    SCAI_LOG_ERROR( logger, "input set: " << *mInputSet )

    SCAI_ASSERT_EQ_ERROR( mInputSet->getGroup(), "LAMAInputSet", "Illegal LAMAInputSet: " << *mInputSet )

    // Now it is safe to cast

    mLAMAInputSet = reinterpret_cast<LAMAInputSet*>( mInputSet.get() );

    const DenseVector<double>& inputX = mLAMAInputSet->getX();
    const CSRSparseMatrix<double>& inputA = mLAMAInputSet->getA();

    // mDistribution = DistributionPtr ( new GenBlockDistribution( inputA.getNumRows(), weight, mComm ) );

    dmemo::DistributionPtr mDistribution;

    mDistribution = inputA.getRowDistributionPtr();

    SCAI_LOG_INFO( logger,
                   "General block distribution of matrix with weight " << weight << ", local size = " << mDistribution->getLocalSize() );

    *mVectorX = inputX;
    mVectorX->redistribute( mDistribution );

    mVectorY->allocate( mDistribution );
    *mVectorY = Scalar( 0 );

    SCAI_LOG_DEBUG( logger, "convert and redistribute the input matrix " << inputA );

    *mMatrixA = inputA;
    mMatrixA->redistribute( mDistribution, mDistribution );

    SCAI_LOG_DEBUG( logger, "matrix " << *mMatrixA << " now available for MV" );

    mMatrixA->setContextPtr( mContext );
    mMatrixA->setCommunicationKind( mCommunicationKind );

    SCAI_LOG_INFO( logger,
                   "Matrix: context at " << *mContext << ", comm = " << mCommunicationKind );
}

void PMVBenchmark::setUp()
{
    mVectorX->prefetch( mContext );
    mMatrixA->prefetch();
    mVectorX->wait();
    mMatrixA->wait();

    SCAI_LOG_INFO( logger,
                   "setUp done for p = " << mComm->getRank() << " : X, A at " << *mContext );
}

void PMVBenchmark::execute()
{
    *mVectorY = *mMatrixA * *mVectorX;
}

void PMVBenchmark::tearDown()
{
    hmemo::ContextPtr host = hmemo::Context::getHostPtr();
    mVectorY->prefetch( host );
    mVectorY->wait();
    SCAI_LOG_INFO( logger, "tearDown done, Y at Host" );
}

void PMVBenchmark::shutdown()
{
    SCAI_ASSERT_ERROR( mLAMAInputSet, "No LAMA input set available" )

    LAMAInputSetComplexityVisitor::getMVComplexity( mLAMAInputSet->getA(), mNumFloatingPointOperations,
                                                    mNumProcessedBytesFloat, mNumProcessedBytesDouble );

    const DenseVector<double>& result = mLAMAInputSet->getY();

    // set the maximal difference that is allowed, depends on ValueType

    Scalar maxDiff = 1E-4f;

    try
    {
        const _DenseVector& computedResult = *mVectorY;

        SCAI_LOG_INFO( logger, "computed result: " << computedResult );

        // Note: there might be type conversion from double to float

        DenseVector<double> correctResult( result );

        SCAI_LOG_INFO( logger, "result ( by inputSet ): " << result );
        SCAI_LOG_INFO( logger, "correct result: " << correctResult );

        Vector& diff = correctResult;

        diff = computedResult - correctResult;

        Scalar diffNorm = maxNorm( diff );

        SCAI_LOG_INFO( logger, "max diff = " << diffNorm );

        SCAI_ASSERT_ERROR( diffNorm < maxDiff, "illegal result" );
    }
    catch ( benchmark::BFException& e )
    {
        std::stringstream message;
        message << e.what() << " std::fabs( result[i] - computedResultValue ) is bigger than " << maxDiff << std::endl;
        throw benchmark::BFException( message.str() );
    }

    mMatrixA.reset();
    mVectorX.reset();
    mVectorY.reset();
}

CounterType PMVBenchmark::getNumFloatingPointOperations() const
{
    return mNumFloatingPointOperations;
}

CounterType PMVBenchmark::getProcessedBytes() const
{
    size_t typeSize = getValueTypeSize();

    if ( typeSize == sizeof( float ) )
    {
        return mNumProcessedBytesFloat;
    }
    else if ( typeSize == sizeof( double ) )
    {
        return mNumProcessedBytesDouble;
    }

    return 0;
}

}

}