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

#include <scai/benchmark/Parser.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

namespace scai
{

using std::istringstream;
using std::map;
using std::string;

namespace lama
{

const string& PMVBenchmark::getGroupId()
{
    static const string id = "PMV";
    return id;
}

PMVBenchmark::PMVBenchmark( const string& argument ) :

    LAMAMPIBenchmark( getCreateId(), getGroupId() ),
    mContext( hmemo::Context::getContextPtr( common::ContextType::Host ) ),
    mNumFloatingPointOperations( 0 ),
    mNumProcessedIndexes( 0 ),
    mNumProcessedValues( 0 )
{
    std::vector<string> argTokens;

    if ( argument == "" )
    {
        tokenize( argTokens, "CSR double", " ,:" );
    }
    else
    {
        tokenize( argTokens, argument, " ,:" );
    }

    SCAI_ASSERT_EQ_ERROR( argTokens.size(), 2, "format, type expected as argument" )

    Format format = str2Format( argTokens[0].c_str() );

    if ( format == Format::UNDEFINED )
    {
        COMMON_THROWEXCEPTION( "unknowm matrix format: " << argTokens[0] )
    }

    mType = common::str2ScalarType( argTokens[1].c_str() );

    if ( mType == common::ScalarType::UNKNOWN )
    {
        COMMON_THROWEXCEPTION( "unknowm value type: " << argTokens[1] )
    }

    // allocate source and target storage of the required type

    mMatrixA.reset( Matrix<ValueType>::getMatrix( format ) );
    mVectorX.reset( new DenseVector<ValueType>() );
    mVectorY.reset( new DenseVector<ValueType>() );

    mArgument = argTokens[0] + ", " + argTokens[1];

    mName = "PMV( " + mArgument + " )";
}

PMVBenchmark::~PMVBenchmark()
{
    // Note: mMatrixA, mVectorX, mVectorY will be freed
}

common::ScalarType PMVBenchmark::getValueType() const
{
    return mType;
}

bool PMVBenchmark::isThreadded() const
{
    return true;
}

const string& PMVBenchmark::getCreateId() const
{
    static string id = createValue();
    return id;
}

const string& PMVBenchmark::getArgument() const
{
    return mArgument;
}

string PMVBenchmark::createValue() 
{
    string value ( "PMatTimesVector" );
    return value;
}

void PMVBenchmark::initialize()
{
    //device initialization + comunicator

    map<string, string> tokens;

    int devNo = -1; // default device

    if ( tokens.count( "CUDA" ) > 0 )
    {
        istringstream tokenStream( tokens["CUDA"] );
        tokenStream >> devNo;
    }

    std::ostringstream idStream;
    idStream << " (Proc " << mComm->getRank() << " LOCAL=";

    mContext = hmemo::Context::getContextPtr();

    //TODO: Aggregate Name at Proc 0 to print it in output

    SCAI_LOG_INFO( logger, "get input set by this id " << mInputSetId );

    // mInputSetId = "inputSetId( argument )" must also be parsed ( "inputSetId", "argument" )

    mInputSet.reset( benchmark::InputSet::createWithArgument( mInputSetId ) );

    SCAI_LOG_ERROR( logger, "input set: " << *mInputSet )

    SCAI_ASSERT_EQ_ERROR( mInputSet->getGroup(), "LAMAInputSet", "Illegal LAMAInputSet: " << *mInputSet )

    // Now it is safe to cast

    mLAMAInputSet = reinterpret_cast<LAMAInputSet*>( mInputSet.get() );

    const DenseVector<double>& inputX = mLAMAInputSet->getX();
    const CSRSparseMatrix<double>& inputA = mLAMAInputSet->getA();

    dmemo::DistributionPtr mDistribution;

    mDistribution = inputA.getRowDistributionPtr();

    *mVectorX = inputX;
    mVectorX->redistribute( mDistribution );

    mVectorY->setSameValue( mDistribution, ValueType( 0 ) );

    SCAI_LOG_DEBUG( logger, "convert and redistribute the input matrix " << inputA );

    *mMatrixA = inputA;
    mMatrixA->redistribute( mDistribution, mDistribution );

    SCAI_LOG_DEBUG( logger, "matrix " << *mMatrixA << " now available for MV" );

    mMatrixA->setContextPtr( mContext );
    mMatrixA->setCommunicationKind( mCommunicationKind );

    SCAI_LOG_ERROR( logger,
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

void PMVBenchmark::getComplexity( 
    CounterType& numFlops,
    CounterType& numProcessedIndexes,
    CounterType& numProcessedValues,
    const Matrix<ValueType>& matrix )
{
    IndexType numRows = matrix.getNumRows(); // used for all formats

    if ( matrix.getMatrixKind() == MatrixKind::DENSE )
    {
        IndexType numCols = matrix.getNumColumns();
        numFlops = numRows * ( 2 * numCols - 1 );
        numProcessedValues = 0;
        numProcessedValues = 2 * numRows * numCols + numRows;
        return;
    }
    else
    {
        IndexType numValues = matrix.getNumValues();
        //every matrix element is multiplied once
        //and all products of one row of the matrix are added together
        //so we have numValues float multiplies and numValues - numRows float adds.
        numFlops = 2 * numValues - numRows;
        //The whole matrix need to be accessed once
        //for each row of the matrix two values of the index array need to be loaded
        numProcessedIndexes = numValues + numRows;
        //we need to load the values of the matrix once
        numProcessedValues  = numValues;
        //we need to load one value of the input vector for each element of
        //the matrix (we ignore the cache)
        numProcessedValues += numValues;
        //we need to write each value of the output vector once
        numProcessedValues += numRows;
        return;
    }
}

void PMVBenchmark::shutdown()
{
    SCAI_ASSERT_ERROR( mLAMAInputSet, "No LAMA input set available" )

    getComplexity( mNumFloatingPointOperations, mNumProcessedIndexes, mNumProcessedValues, mLAMAInputSet->getA() );

    const DenseVector<double>& result = mLAMAInputSet->getY();

    // set the maximal difference that is allowed, depends on ValueType

    RealType<ValueType> maxDiff = 1E-4f;

    try
    {
        const DenseVector<ValueType>& computedResult = *mVectorY;

        SCAI_LOG_INFO( logger, "computed result: " << computedResult );

        DenseVector<ValueType> correctResult( result );

        SCAI_LOG_INFO( logger, "result ( by inputSet ): " << result );
        SCAI_LOG_INFO( logger, "correct result: " << correctResult );

        Vector<ValueType>& diff = correctResult;

        diff = computedResult - correctResult;

        auto diffNorm = maxNorm( diff );

        SCAI_LOG_INFO( logger, "max diff = " << diffNorm );

        SCAI_ASSERT_ERROR( diffNorm < maxDiff, "illegal result" );
    }
    catch ( benchmark::BFException& e )
    {
        std::ostringstream message;
        message << e.what() << " std::fabs( result[i] - computedResultValue ) is bigger than " << maxDiff << std::endl;
        throw benchmark::BFException( message.str() );
    }

    mMatrixA.reset();
    mVectorX.reset();
    mVectorY.reset();
}

CounterType PMVBenchmark::getNumFloatingPointOperations() const
{
    // ToDo: for complex values we have some more floating point operations

    return mNumFloatingPointOperations;
}

CounterType PMVBenchmark::getProcessedBytes() const
{
    return   mNumProcessedValues * common::typeSize( getValueType() )
           + mNumProcessedIndexes * sizeof( IndexType );
}

}

}
