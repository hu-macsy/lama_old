/**
 * @file Matrix.cpp
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
 * @brief Matrix.cpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */

// hpp
#include <lama/matrix/Matrix.hpp>

// others
#include <lama/DenseVector.hpp>
#include <lama/distribution/NoDistribution.hpp>

#include <lama/exception/LAMAAssert.hpp>

using namespace boost;

namespace lama
{

LAMA_LOG_DEF_LOGGER( Matrix::logger, "Matrix" );

Matrix::Matrix( const Matrix& other )

    : Distributed( other ), mColDistribution( other.mColDistribution ), mNumRows( other.mNumRows ), mNumColumns(
        other.mNumColumns ), mCommunicationKind( other.mCommunicationKind )
{
    LAMA_LOG_INFO( logger, "Creating copy of " << other << " with same distributions." );
}

Matrix::Matrix( const Matrix& other, DistributionPtr distribution, DistributionPtr colDistribution )

    : Distributed( distribution ), mColDistribution( colDistribution ), mNumRows( other.mNumRows ), mNumColumns(
        other.mNumColumns ), mCommunicationKind( other.mCommunicationKind )
{
    // Very important: here we check that new distributions fit the matrix

    checkSettings();

    LAMA_LOG_INFO( logger,
                   "Creating copy of " << other << " with new distributions: " << "row = " << getDistribution() << ", col = " << getColDistribution() );

}

Matrix::Matrix( const IndexType numRows, const IndexType numColumns )

    : Distributed( DistributionPtr( new NoDistribution( numRows ) ) ), mColDistribution(
        DistributionPtr( new NoDistribution( numColumns ) ) ), mNumRows( numRows ), mNumColumns(
            numColumns )
{
    setDefaultKind();

    LAMA_LOG_INFO( logger, "Creating a replicated Matrix of size " << mNumRows << " x " << mNumColumns );
}

/* ----------------------------------------------------------------------- */

void Matrix::checkSettings() const
{
    if( !mColDistribution )
    {
        LAMA_THROWEXCEPTION( "NULL pointer for column distribution" );
    }

    if( mNumRows != getDistribution().getGlobalSize() )
    {
        LAMA_THROWEXCEPTION(
            "row distribution " << getDistribution() << ": global size mismatches #rows = " << mNumRows );
    }

    if( mNumColumns != getColDistribution().getGlobalSize() )
    {
        LAMA_THROWEXCEPTION(
            "col distribution " << getColDistribution() << ": global size mismatches #columns = " << mNumColumns );
    }
}

Matrix::Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )

    : Distributed( rowDistribution )
{
    setDistributedMatrix( rowDistribution, colDistribution );

    setDefaultKind();

    LAMA_LOG_INFO( logger,
                   "Construct a Matrix of size " << mNumRows << " x " << mNumColumns << " with the distribution " << getDistribution() );

    checkSettings();
}

Matrix::Matrix( DistributionPtr distribution )
    : Distributed( distribution )
{
    setDistributedMatrix( distribution, distribution );

    LAMA_LOG_INFO( logger,
                   "Construct a square Matrix of size " << mNumRows << " x " << mNumColumns << " with the row/col distribution " << getDistribution() );

    checkSettings();
}

Matrix::Matrix()

    : Distributed( DistributionPtr( new NoDistribution( 0 ) ) ), mColDistribution(
        DistributionPtr( new NoDistribution( 0 ) ) ), mNumRows( 0 ), mNumColumns( 0 )
{
    setDefaultKind();
}

Matrix::~Matrix()
{
    LAMA_LOG_DEBUG( logger, "~Matrix" );
}

void Matrix::setDefaultKind()
{
    mCommunicationKind = ASYNCHRONOUS;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setDistributedMatrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_ASSERT_ERROR( rowDistribution, "NULL row distribution for matrix not allowd" );
    LAMA_ASSERT_ERROR( colDistribution, "NULL column distribution for matrix not allowd" );
    mNumRows = rowDistribution->getGlobalSize();
    setDistributionPtr( rowDistribution );
    mNumColumns = colDistribution->getGlobalSize();
    mColDistribution = colDistribution;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setReplicatedMatrix( const IndexType numRows, const IndexType numColumns )
{
    DistributionPtr rowDist( new NoDistribution( numRows ) );

    if( numRows == numColumns )
    {
        setDistributedMatrix( rowDist, rowDist );
    }
    else
    {
        setDistributedMatrix( rowDist, DistributionPtr( new NoDistribution( numColumns ) ) );
    }
}

/* ---------------------------------------------------------------------------------*/

std::auto_ptr<_LAMAArray> Matrix::createArray() const
{
    // Static method of _LAMAArray provides exactly the needed functionality.

    return _LAMAArray::create( getValueType() );
}

/* ---------------------------------------------------------------------------------*/

std::auto_ptr<Matrix> Matrix::create( const IndexType numRows, const IndexType numColumns ) const
{
    std::auto_ptr<Matrix> matrix = create();
    matrix->allocate( numRows, numColumns );
    return matrix;
}

/* ---------------------------------------------------------------------------------*/

std::auto_ptr<Matrix> Matrix::create( const IndexType size ) const
{
    std::auto_ptr<Matrix> matrix = create();
    DistributionPtr dist( new NoDistribution( size ) );
    matrix->allocate( dist, dist );
    return matrix;
}

/* ---------------------------------------------------------------------------------*/

std::auto_ptr<Matrix> Matrix::create( DistributionPtr rowDistribution, DistributionPtr colDistribution ) const
{
    std::auto_ptr<Matrix> matrix = create();
    matrix->allocate( rowDistribution, colDistribution );
    return matrix;
}

/* ---------------------------------------------------------------------------------*/

std::auto_ptr<Matrix> Matrix::create( DistributionPtr distribution ) const
{
    std::auto_ptr<Matrix> matrix = create();
    matrix->allocate( distribution, distribution );
    return matrix;
}

/* ---------------------------------------------------------------------------------*/

VectorPtr Matrix::createDenseVector( DistributionPtr distribution, const Scalar value ) const
{
    Scalar::ScalarType matrixValueType = getValueType();

    LAMA_LOG_INFO( logger, "create vector of type " << matrixValueType );

    switch( matrixValueType )
    {
    case Scalar::DOUBLE:
        return VectorPtr( new DenseVector<double>( distribution, value.getValue<double>() ) );
    case Scalar::FLOAT:
        return VectorPtr( new DenseVector<float>( distribution, value.getValue<float>() ) );
    default:
        LAMA_THROWEXCEPTION( "unsupported vector type : " << matrixValueType );
    }
}

/* ---------------------------------------------------------------------------------*/

Scalar Matrix::operator()( IndexType i, IndexType j ) const
{
    return getValue( i, j );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setCommunicationKind( SyncKind communicationKind )
{
    mCommunicationKind = communicationKind;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setContext( ContextPtr localContext, ContextPtr haloContext )
{
    // default implementation for matrices that do not support halo context

    if( *localContext != *haloContext )
    {
        LAMA_LOG_WARN( logger, *this << ": halo context = " << *haloContext << " ignored" );
    }

    setContext( localContext );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::inheritAttributes( const Matrix& other )
{
    setCommunicationKind( other.getCommunicationKind() );
    setContext( other.getContextPtr() );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::writeAt( std::ostream& stream ) const
{
    stream << "Matrix(" << mNumRows << "x" << mNumColumns << ")";
}

/* ---------------------------------------------------------------------------------*/

void Matrix::matrix2CSRGraph(
    IndexType* /*xadj*/,
    IndexType* /*adjncy*/,
    IndexType* /*vwgt*/,
    CommunicatorPtr /*comm*/,
    const IndexType* /*globalRowIndices*/,
    IndexType* /*vtxdist = NULL*/) const
{
    LAMA_THROWEXCEPTION( "Transformation of "<< *this << "to CSR Graph not specialized." );
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator=( const Matrix& other )
{
    // assignment operator is just implemented by the assign method

    LAMA_LOG_INFO( logger, *this << ": operator = " << other );

    this->assign( other );

    LAMA_LOG_INFO( logger, *this << ": end operator = " << other );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

/**
 * @brief the assignment operator for a scalar matrix multiplication.
 */
Matrix& Matrix::operator=( Expression<Scalar,Matrix,Times> exp )
{
    //TODO Implement --> SCAL --> VectorInterface?
    //LAMA_THROWEXCEPTION("Assignement operator for Scalar * Matrix is not yet implemented")

    const Matrix& A = exp.getArg2();
    const Scalar& s = exp.getArg1();
    this->matrixTimesScalar( A, s );
    return *this;
}

/**
 * @brief the assignment operator for a matrix matrix multiplication.
 */
Matrix& Matrix::operator=( Expression<Matrix,Matrix,Times> exp )
{
    Scalar zero( 0 );
    Scalar one( 1 );
    Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> exp1( one, exp );
    Expression<Scalar,Matrix,Times> exp2( zero, exp.getArg1() );
    *this = Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
                exp1, exp2 );
    return *this;
}

/**
 * @brief the assignment operator for a matrix matrix multiplication.
 */
Matrix& Matrix::operator=( Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> exp )
{
    Scalar zero( 0 );
    Expression<Scalar,Matrix,Times> exp1( zero, exp.getArg2().getArg1() );
    *this = Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
                exp, exp1 );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::swapMatrix( Matrix& other )
{
    Distributed::swap( other );
    std::swap( mNumRows, other.mNumRows );
    std::swap( mNumColumns, other.mNumColumns );
    std::swap( mColDistribution, other.mColDistribution );
}

/* ---------------------------------------------------------------------------------*/

double Matrix::getSparsityRate() const
{
    return (double) getNumValues() / mNumRows / mNumColumns;
}

/* ---------------------------------------------------------------------------------*/

/**
 * @brief the assignment operator for a GEMM expression.
 */
Matrix& Matrix::operator=(
    Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    LAMA_LOG_INFO( logger, "operator=:  A * B * alpha + C * beta " );

    const Matrix& A = exp.getArg1().getArg2().getArg1();
    const Matrix& B = exp.getArg1().getArg2().getArg2();
    const Matrix& C = exp.getArg2().getArg2();
    const Scalar& alpha = exp.getArg1().getArg1();
    const Scalar& beta = exp.getArg2().getArg1();
    const Scalar zero( 0.0 );

    //Do sanity checks

    const Distribution& rowDistA = A.getDistribution();
    const Distribution& colDistA = A.getColDistribution();
    const Distribution& rowDistB = B.getDistribution();
    const Distribution& colDistB = B.getColDistribution();
    const Distribution& rowDistC = C.getDistribution();
    const Distribution& colDistC = C.getColDistribution();

    if( colDistA != rowDistB )
    {
        LAMA_THROWEXCEPTION(
            "Distribution of " << A << " = " << A.getColDistribution() << " does not match distribution of " << B << " = " << B.getDistribution() );
    }

    if( rowDistA != rowDistC && beta != zero )
    {
        LAMA_THROWEXCEPTION(
            "Distribution of " << A << " = " << A.getDistribution() << " does not match distribution of " << C << " = " << C.getDistribution() );
    }

    if( colDistB != colDistC && beta != zero )
    {
        LAMA_THROWEXCEPTION(
            "Distribution of " << B << " = " << B.getColDistribution() << " does not match distribution of " << C << " = " << C.getColDistribution() );
    }

    // lhs matrix will be allocated with ( rowDistA, colDistB )

    //size checks are needed because NoDistribution does not compare sizes

    if( A.getNumColumns() != B.getNumRows() )
    {
        LAMA_THROWEXCEPTION(
            "Number of rows of " << A << " = " << A.getNumRows() << " does not match the number of columns of " << B << " = " << B.getNumColumns() );
    }
    if( A.getNumRows() != C.getNumRows() && beta != zero )
    {
        LAMA_THROWEXCEPTION(
            "Number of rows of " << A << " = " << A.getNumRows() << " does not match the number of rows of " << C << " = " << C.getNumRows() );
    }
    if( B.getNumColumns() != C.getNumColumns() && beta != zero )
    {
        LAMA_THROWEXCEPTION(
            "Number of columns of " << B << " = " << B.getNumColumns() << " does not match the number of columns of " << C << " = " << C.getNumColumns() );
    }

    LAMA_LOG_INFO( logger, "Context of this before matrixTimesMatrix = " << this->getContext() );

    A.matrixTimesMatrix( *this, alpha, B, beta, C );

    LAMA_LOG_INFO( logger, "end operator=:  A * B * alpha + C * beta " );

    LAMA_LOG_INFO( logger, "Context of this after matrixTimesMatrix = " << this->getContext() );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

/**
 * @brief the assignment operator for a MM addition.
 */
Matrix& Matrix::operator=( Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    LAMA_LOG_INFO( logger, "operator=:  A * alpha + B * beta " );

    const Matrix& A = exp.getArg1().getArg2();
    const Matrix& B = exp.getArg2().getArg2();
    const Scalar& alpha = exp.getArg1().getArg1();
    const Scalar& beta = exp.getArg2().getArg1();
    const Scalar zero( 0.0 );

    if( beta == zero )
    {
        this->matrixTimesScalar( A, alpha );
        return *this;
    }

    if( alpha == zero )
    {
        this->matrixTimesScalar( B, beta );
        return *this;
    }

    // Do sanity checks

    const Distribution& rowDistA = A.getDistribution();
    const Distribution& colDistA = A.getColDistribution();
    const Distribution& rowDistB = B.getDistribution();
    const Distribution& colDistB = B.getColDistribution();

    if( rowDistA != rowDistB )
    {
        LAMA_THROWEXCEPTION(
            "Row distribution of " << A << " = " << rowDistA << " does not match distribution of " << B << " = " << rowDistB );
    }

    if( colDistA != colDistB )
    {
        LAMA_THROWEXCEPTION(
            "Column distribution of " << A << " = " << colDistA << " does not match distribution of " << B << " = " << colDistB );
    }

    this->matrixPlusMatrix( alpha, A, beta, B );

    return *this;
}

}
