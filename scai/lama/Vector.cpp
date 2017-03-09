/**
 * @file Vector.cpp
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
 * @brief Implementations of methods for class Vector.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/Vector.hpp>

// local library

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Distribution.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/io/PartitionIO.hpp>

// tracing
#include <scai/tracing.hpp>

// std
#include <map>
#include <ostream>

using namespace scai::common;
using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( Vector::logger, "Vector" )

/* ---------------------------------------------------------------------------------- */

const char* _Vector::kind2Str( const VectorKind vectorKind )
{
    switch ( vectorKind )
    {
        case DENSE:
            return "DENSE";
            break;

        case SPARSE:
            return "SPARSE";
            break;

        case UNDEFINED:
            return "Undefined";
            break;
    }

    return "<illegal_vector_kind>";
}

_Vector::VectorKind _Vector::str2Kind( const char* str )

{
    for ( int kind = DENSE; kind < UNDEFINED; ++kind )
    {
        if ( strcmp( kind2Str( VectorKind( kind ) ), str ) == 0 )
        {
            return VectorKind( kind );
        }
    }

    return UNDEFINED;
}

/* ---------------------------------------------------------------------------------------*/
/*    VectorKind opertor<<                                                                */
/* ---------------------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const _Vector::VectorKind& kind )
{
    stream << _Vector::kind2Str( kind );
    return stream;
}

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a vector                                                          */
/* ---------------------------------------------------------------------------------------*/

Vector* Vector::getVector( const VectorKind kind, const common::scalar::ScalarType valueType )
{
    VectorCreateKeyType vectype( kind, valueType );
    return Vector::create( vectype );
}

Vector* Vector::getDenseVector(
    const common::scalar::ScalarType valueType,
    DistributionPtr distribution,
    ContextPtr context )
{
    VectorCreateKeyType vectype( Vector::DENSE, valueType );
    Vector* v = Vector::create( vectype );
    v->allocate( distribution );

    if ( context )
    {
        v->setContextPtr( context );
    }

    return v;
}

/* ---------------------------------------------------------------------------------------*/
/*    Constructor / Destructor                                                            */
/* ---------------------------------------------------------------------------------------*/

Vector::Vector( const IndexType size, hmemo::ContextPtr context ) :

    Distributed( DistributionPtr( new NoDistribution( size ) ) ),
    mContext( context )
{
    if ( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger, "Vector(" << size << "), replicated, on " << *mContext )
}

Vector::Vector( DistributionPtr distribution, hmemo::ContextPtr context )
    : Distributed( distribution ), mContext( context )
{
    if ( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger,
                   "Vector(" << distribution->getGlobalSize() << ") with " << getDistribution() << " constructed" )
}

Vector::Vector( const Vector& other )
    : Distributed( other ), mContext( other.getContextPtr() )
{
    SCAI_ASSERT_ERROR( mContext, "NULL context not allowed" )
    SCAI_LOG_INFO( logger, "Vector(" << other.getDistribution().getGlobalSize() << "), distributed, copied" )
}

Vector::~Vector()
{
    SCAI_LOG_DEBUG( logger, "~Vector(" << getDistribution().getGlobalSize() << ")" )
}

/* ---------------------------------------------------------------------------------------*/
/*    Reading vector from a file, only host reads                                         */
/* ---------------------------------------------------------------------------------------*/

void Vector::readFromSingleFile( const std::string& fileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const PartitionId MASTER = 0;
    const PartitionId myRank = comm->getRank();

    // this is a bit tricky stuff, but it avoids an additional copy from array -> vector

    IndexType globalSize;

    if ( myRank == MASTER )
    {
        SCAI_LOG_INFO( logger, *comm << ": read local array from file " << fileName )

        globalSize = readLocalFromFile( fileName );
    }

    comm->bcast( &globalSize, 1, MASTER );

    if ( myRank != MASTER )
    {
        clearValues();
    }

    DistributionPtr distribution( new SingleDistribution( globalSize, comm, MASTER ) );

    setDistributionPtr( distribution );

    SCAI_LOG_INFO( logger, "readFromSingleFile, vector = " << *this )
}

/* ---------------------------------------------------------------------------------------*/
/*    Reading vector from a file, every processor reads its partition                     */
/* ---------------------------------------------------------------------------------------*/

void Vector::readFromSingleFile( const std::string& fileName, const DistributionPtr distribution )
{
    if ( distribution.get() == NULL )
    {
        SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << ", master only" )
        readFromSingleFile( fileName );
        return;
    }

    const IndexType n = distribution->getBlockDistributionSize();

    if ( n == nIndex )
    {
        SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << " ), master only + redistribute" )
        readFromSingleFile( fileName );
        redistribute( distribution );
        return;
    }

    // we have a block distribution, so every processor reads its own part

    IndexType first = 0;

    if ( n > 0 )
    {
        first = distribution->local2global( 0 );   // first global index
    }

    bool error = false;

    SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << " ), block dist = " << *distribution 
                           << ", read my block, first = " << first << ", n = " << n )

    try
    {
        IndexType localSize = readLocalFromFile( fileName, first, n );
        error = localSize != distribution->getLocalSize();
    }
    catch ( Exception& ex )
    {
        SCAI_LOG_ERROR( logger, ex.what() )
        error = true;
    }

    error = distribution->getCommunicator().any( error );

    if ( error )
    {
        COMMON_THROWEXCEPTION( "readFromSingleFile " << fileName << " failed, dist = " << *distribution )
    }

    setDistributionPtr( distribution );
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromPartitionedFile( const std::string& myPartitionFileName, DistributionPtr dist )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    bool errorFlag = false;

    IndexType localSize = 0;

    try
    {
        localSize = readLocalFromFile( myPartitionFileName );

        if ( dist.get() )
        {
            // size of storage must match the local size of distribution

            SCAI_ASSERT_EQUAL( dist->getLocalSize(), localSize, "serious mismatch: local matrix in " << myPartitionFileName << " has illegal local size" )
        }
    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_ERROR( logger, *comm << ": failed to read " << myPartitionFileName << ": " << e.what() )
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "error reading partitioned matrix" )
    }

    DistributionPtr vectorDist = dist;

    if ( !vectorDist.get() )
    {
        // we have no distribution so assume a general block distribution

        IndexType globalSize = comm->sum( localSize );

        vectorDist.reset( new GenBlockDistribution( globalSize, localSize, comm ) );
    }

    setDistributionPtr( vectorDist );   // distribution matches size of local part
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromFile( const std::string& vectorFileName, const std::string& distributionFileName )
{
    // read the distribution

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr distribution;

    if ( distributionFileName == "BLOCK" )
    {
        // for a single file we set a BlockDistribution

        if ( vectorFileName.find( "%r" ) == std::string::npos )
        {
            PartitionId root = 0;

            IndexType numRows = nIndex;

            if ( comm->getRank() == root )
            {
                numRows = FileIO::getStorageSize( vectorFileName );
            }

            comm->bcast( &numRows, 1, root );

            distribution.reset( new BlockDistribution( numRows, comm ) );
        }

        // for a partitioned file general block distribution is default
    }
    else
    {
        distribution = PartitionIO::readDistribution( distributionFileName, comm );
    }

    readFromFile( vectorFileName, distribution );
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromFile( const std::string& fileName, DistributionPtr distribution )
{
    SCAI_LOG_INFO( logger, *this << ": readFromFile( " << fileName << " )" )

    std::string newFileName = fileName;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // take default

    if ( distribution.get() )
    {
        comm = distribution->getCommunicatorPtr();
    }

    bool isPartitioned;

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, *comm );

    if ( !isPartitioned )
    {
        readFromSingleFile( newFileName, distribution );
    }
    else
    {
        readFromPartitionedFile( newFileName, distribution );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*    Assignment operator                                                                 */
/* ---------------------------------------------------------------------------------------*/

Vector& Vector::operator=( const Expression_SV_SV& expression )
{
    SCAI_LOG_DEBUG( logger, "this = a * vector1 + b * vector2, check vector1.size() == vector2.size()" )

    const Scalar& alpha = expression.getArg1().getArg1();
    const Scalar& beta  = expression.getArg2().getArg1();
    const Vector& x = expression.getArg1().getArg2();
    const Vector& y = expression.getArg2().getArg2();

    // Note: all checks are done the vector specific implementations

    vectorPlusVector( alpha, x, beta, y );

    return *this;
}

Vector& Vector::operator=( const Expression_SMV& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * matrix * vectorX -> this = alpha * matrix * vectorX + 0.0 * this" )
    const Scalar beta( 0.0 );
    Expression_SV exp2( beta, *this );
    Expression_SMV_SV tmpExp( expression, exp2 );
    const Vector& vectorX = expression.getArg2().getArg2();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix& matrix = expression.getArg2().getArg1();
        DistributionPtr dist = matrix.getRowDistributionPtr();
        allocate( dist );
        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

Vector& Vector::operator=( const Expression_SVM& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * vectorX * matrix -> this = alpha * vectorX * matrix + 0.0 * this" )
    const Scalar beta( 0.0 );
    Expression_SV exp2( beta, *this );
    Expression_SVM_SV tmpExp( expression, exp2 );
    const Vector& vectorX = expression.getArg2().getArg1();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix& matrix = expression.getArg2().getArg2();
        DistributionPtr dist = matrix.getColDistributionPtr();
        allocate( dist );
        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

Vector& Vector::operator=( const Expression_SMV_SV& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SMV_SV )" )
    const Expression_SMV& exp1 = expression.getArg1();
    const Expression_SV& exp2 = expression.getArg2();
    const Scalar& alpha = exp1.getArg1();
    const Expression<Matrix, Vector, Times>& matrixTimesVectorExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();
    const Matrix& matrix = matrixTimesVectorExp.getArg1();
    const Vector& vectorX = matrixTimesVectorExp.getArg2();
    Vector* resultPtr = this;
    VectorPtr tmpResult;

    if ( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult.reset( Vector::create( this->getCreateValue() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call matrixTimesVector with matrix = " << matrix )
    matrix.matrixTimesVector( *resultPtr, alpha, vectorX, beta, vectorY );

    if ( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression_SVM_SV& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SVM_SV )" )
    const Expression_SVM& exp1 = expression.getArg1();
    const Expression_SV& exp2 = expression.getArg2();
    const Scalar& alpha = exp1.getArg1();
    const Expression<Vector, Matrix, Times>& vectorTimesMatrixExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();
    const Vector& vectorX = vectorTimesMatrixExp.getArg1();
    const Matrix& matrix = vectorTimesMatrixExp.getArg2();
    Vector* resultPtr = this;
    VectorPtr tmpResult;

    if ( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult.reset( Vector::create( this->getCreateValue() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call vectorTimesMatrix with matrix = " << matrix )
    matrix.vectorTimesMatrix( *resultPtr, alpha, vectorX, beta, vectorY );

    if ( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression_SV& expression )
{
    const Scalar& alpha = expression.getArg1();
    const Vector& x = expression.getArg2();

    vectorPlusVector( alpha, x, 0, x );

    return *this;
}

Vector& Vector::operator=( const Expression_VV& expression )
{
    SCAI_LOG_DEBUG( logger, "operator=, SVV( alpha, x, y) -> x * y" )

    const Vector& x = expression.getArg1();
    const Vector& y = expression.getArg2();

    Scalar alpha( 1 );

    vectorTimesVector( alpha, x, y );

    return *this;
}

Vector& Vector::operator=( const Expression_SVV& expression )
{
    const Scalar& alpha = expression.getArg1();

    const Expression_VV& exp = expression.getArg2();
    const Vector& x = exp.getArg1();
    const Vector& y = exp.getArg2();

    vectorTimesVector( alpha, x, y );

    return *this;
}

Vector& Vector::operator=( const Vector& other )
{
    assign( other );

    return *this;
}

Vector& Vector::operator=( const Expression_SV_S& expression )
{
    const Expression_SV& exp = expression.getArg1();
    const Scalar& alpha = exp.getArg1();
    const Vector& x = exp.getArg2();
    const Scalar& beta = expression.getArg2();

    vectorPlusScalar( alpha, x, beta );

    return *this;
}

Vector& Vector::operator=( const Scalar value )
{
    assign( value );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/
/*   Compound assignments                                                                 */
/* ---------------------------------------------------------------------------------------*/

Vector& Vector::operator*=( const Scalar value )
{
    return operator=( Expression_SV( value, *this ) );
}

Vector& Vector::operator*=( const Vector& other )
{
    scale( other );
    return *this;
}

Vector& Vector::operator/=( const Scalar value )
{
    Expression<Scalar, Vector, Times> exp1( Scalar( 1.0 ) / value, *this );
    return operator=( exp1 );
}

Vector& Vector::operator+=( const Vector& other )
{
    return operator=( Expression_SV_SV( Expression_SV( Scalar( 1 ), other ), Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Vector& other )
{
    return operator=( Expression_SV_SV( Expression_SV( Scalar( 1 ), *this ), Expression_SV( Scalar( -1 ), other ) ) );
}

Vector& Vector::operator+=( const Scalar value )
{
    add( value );
    return *this;
}

Vector& Vector::operator-=( const Scalar value )
{
    add( -value );
    return *this;
}

Vector& Vector::operator+=( const Expression_SV& exp )
{
    return operator=( Expression_SV_SV( exp, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Expression_SV& exp )
{
    Expression_SV minusExp( -exp.getArg1(), exp.getArg2() );
    return operator=( Expression_SV_SV( minusExp, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Expression_SMV& expression )
{
    return operator=( Expression_SMV_SV( expression, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Expression_SVM& expression )
{
    return operator=( Expression_SVM_SV( expression, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Expression_SMV& exp )
{
    Expression_SMV minusExp( -exp.getArg1(), exp.getArg2() );
    return operator=( Expression_SMV_SV( minusExp, Expression_SV( Scalar( 1 ), *this ) ) );
}

/* ---------------------------------------------------------------------------------------*/
/*   assign operations                                                                    */
/* ---------------------------------------------------------------------------------------*/

void Vector::assign( const Vector& other )
{
    SCAI_LOG_INFO( logger, "assign other = " << other )

    setDistributionPtr( other.getDistributionPtr() );

    switch ( other.getVectorKind() )
    {
        case Vector::DENSE:
        {
            const _DenseVector& denseOther = reinterpret_cast<const _DenseVector&>( other );
            setDenseValues( denseOther.getLocalValues() );
            break;
        }
        case Vector::SPARSE:
        {
            const _SparseVector& sparseOther = reinterpret_cast<const _SparseVector&>( other );
            setSparseValues( sparseOther.getNonZeroIndexes(), sparseOther.getNonZeroValues(), sparseOther.getZero() );
            break;
        }
        default:

            COMMON_THROWEXCEPTION( "illegal vector kind, other = " << other.getVectorKind() )
    }
}

void Vector::assign( const _HArray& localValues, DistributionPtr dist )
{
    SCAI_ASSERT_EQ_ERROR( localValues.size(), dist->getLocalSize(), "Mismatch local size of vecotr" )

    setDistributionPtr( dist );
    setDenseValues( localValues );
}

void Vector::assign( const _HArray& globalValues )
{
    SCAI_LOG_INFO( logger, "assign vector with globalValues = " << globalValues )

    setDistributionPtr( DistributionPtr( new NoDistribution( globalValues.size() ) ) );
    setDenseValues( globalValues );
}

/* ---------------------------------------------------------------------------------------*/
/*   writeToFile                                                                          */
/* ---------------------------------------------------------------------------------------*/

void Vector::writeToSingleFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType,
    const FileIO::FileMode fileMode
) const
{
    if ( getDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        if ( comm->getRank() == 0 )
        {
            writeLocalToFile( fileName, fileType, dataType, fileMode );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file

        comm->synchronize();
    }
    else
    {
        // writing a distributed vector into a single file requires redistributon

        common::unique_ptr<Vector> repV( copy() );
        repV->replicate();
        repV->writeToSingleFile( fileName, fileType, dataType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------------*/

void Vector::replicate()
{
    if ( getDistribution().isReplicated() )
    {
        return;
    }

    redistribute( DistributionPtr( new NoDistribution( size() ) ) );
}

/* ---------------------------------------------------------------------------------------*/

void Vector::writeToPartitionedFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType,
    const FileIO::FileMode fileMode ) const
{
    bool errorFlag = false;

    try
    {
        writeLocalToFile( fileName, fileType, dataType, fileMode );
    }
    catch ( common::Exception& e )
    {
        errorFlag = true;
    }

    const Communicator& comm = getDistribution().getCommunicator();

    errorFlag = comm.any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Partitioned IO of vector failed" )
    }
}

/* ---------------------------------------------------------------------------------------*/

void Vector::writeToFile(
    const std::string& fileName,
    const std::string& fileType,               /* = "", take IO type by suffix   */
    const common::scalar::ScalarType dataType, /* = UNKNOWN, take defaults of IO type */
    const FileIO::FileMode fileMode            /* = DEFAULT_MODE */
) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    std::string newFileName = fileName;

    bool writePartitions;

    const Communicator& comm = getDistribution().getCommunicator();

    PartitionIO::getPartitionFileName( newFileName, writePartitions, comm );

    if ( !writePartitions )
    {
        writeToSingleFile( newFileName, fileType, dataType, fileMode );
    }
    else
    {
        // matrix_%r.mtx -> matrix_0.4.mtx,  ..., matrix_3.4.mtxt
        writeToPartitionedFile( newFileName, fileType, dataType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*   Miscellaneous                                                                        */
/* ---------------------------------------------------------------------------------------*/

void Vector::swapVector( Vector& other )
{
    // swaps only on this base class, not whole vectors
    mContext.swap( other.mContext );
    Distributed::swap( other );
}

void Vector::writeAt( std::ostream& stream ) const
{
    stream << "Vector(" << getDistributionPtr()->getGlobalSize() << ")";
}

void Vector::setContextPtr( ContextPtr context )
{
    SCAI_ASSERT_DEBUG( context, "NULL context invalid" )

    if ( mContext->getType() != context->getType() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new context = " << *context << ", old context = " << *mContext )
    }

    mContext = context;
}

void Vector::prefetch() const
{
    prefetch( mContext );
}

/* ---------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
