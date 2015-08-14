/**
 * @file Matrix.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @since 1.0.0
 */

// hpp
#include <scai/lama/matrix/Matrix.hpp>

// others
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>

#include <scai/lama/exception/LAMAAssert.hpp>

using namespace common;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( Matrix::logger, "Matrix" )

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a matrix                                                          */
/* ---------------------------------------------------------------------------------------*/

Matrix* Matrix::getMatrix( const MatrixStorageFormat format, common::ScalarType type )
{
    using ::operator<<;   // important to make operator<< outside this namespace available

    MatrixCreateKeyType val( format, type );
    SCAI_LOG_INFO( logger, "getMatrix uses Factory::create " << val )
    return create( val );
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( const Matrix& other ) : 

    Distributed( other ), 
    mColDistribution( other.mColDistribution ), 
    mNumRows( other.mNumRows ), 
    mNumColumns( other.mNumColumns ), 
    mCommunicationKind( other.mCommunicationKind )
{
    SCAI_LOG_INFO( logger, "Creating copy of " << other << " with same distributions." )
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist ) :

    Distributed( rowDist ),
    mColDistribution( colDist ),
    mNumRows( other.mNumRows ),
    mNumColumns( other.mNumColumns ),
    mCommunicationKind( other.mCommunicationKind )
{
    // Very important: here we check that new distributions fit the matrix

    checkSettings();

    SCAI_LOG_INFO( logger,
                   "Creating copy of " << other << " with new distributions: " << "row = " << getDistribution() << ", col = " << getColDistribution() )
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( const IndexType numRows, const IndexType numColumns ) :

    Distributed( DistributionPtr( new NoDistribution( numRows ) ) ),
    mColDistribution( DistributionPtr( new NoDistribution( numColumns ) ) ),
    mNumRows( numRows ),
    mNumColumns( numColumns )
{
    setDefaultKind();

    SCAI_LOG_INFO( logger, "Creating a replicated Matrix of size " << mNumRows << " x " << mNumColumns )
}

/* ----------------------------------------------------------------------- */

void Matrix::setIdentity( const IndexType n )
{
    // take replicated distribution and use pure method

    setIdentity( DistributionPtr( new NoDistribution( n ) ) );
}

/* ----------------------------------------------------------------------- */

void Matrix::checkSettings() const
{
    if ( !mColDistribution )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for column distribution" )
    }

    if ( mNumRows != getDistribution().getGlobalSize() )
    {
        COMMON_THROWEXCEPTION(
            "row distribution " << getDistribution() << ": global size mismatches #rows = " << mNumRows );
    }

    if ( mNumColumns != getColDistribution().getGlobalSize() )
    {
        COMMON_THROWEXCEPTION(
            "col distribution " << getColDistribution() << ": global size mismatches #columns = " << mNumColumns );
    }
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
    : Distributed( rowDistribution )
{
    setDistributedMatrix( rowDistribution, colDistribution );

    setDefaultKind();

    SCAI_LOG_INFO( logger,
                   "Construct a Matrix of size " << mNumRows << " x " << mNumColumns << " with the distribution " << getDistribution() )

    checkSettings();
}

Matrix::Matrix( DistributionPtr distribution )
    : Distributed( distribution )
{
    setDistributedMatrix( distribution, distribution );

    SCAI_LOG_INFO( logger,
                   "Construct a square Matrix of size " << mNumRows << " x " << mNumColumns << " with the row/col distribution " << getDistribution() )

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
    SCAI_LOG_DEBUG( logger, "~Matrix" )
}

void Matrix::setDefaultKind()
{
    mCommunicationKind = ASYNCHRONOUS;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setDistributedMatrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_ASSERT_ERROR( rowDistribution, "NULL row distribution for matrix not allowed" )
    SCAI_ASSERT_ERROR( colDistribution, "NULL column distribution for matrix not allowed" )
    mNumRows = rowDistribution->getGlobalSize();
    setDistributionPtr( rowDistribution );
    mNumColumns = colDistribution->getGlobalSize();
    mColDistribution = colDistribution;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setReplicatedMatrix( const IndexType numRows, const IndexType numColumns )
{
    DistributionPtr rowDist( new NoDistribution( numRows ) );

    if ( numRows == numColumns )
    {
        setDistributedMatrix( rowDist, rowDist );
    }
    else
    {
        setDistributedMatrix( rowDist, DistributionPtr( new NoDistribution( numColumns ) ) );
    }
}

/* ---------------------------------------------------------------------------------*/

ContextArray* Matrix::createArray() const
{
    // Static method of ContextArray provides exactly the needed functionality.

    return ContextArray::create( getValueType() );
}

/* ---------------------------------------------------------------------------------*/

Vector* Matrix::createDenseVector( DistributionPtr distribution, const Scalar value ) const
{
    common::ScalarType matrixValueType = getValueType();

    SCAI_LOG_INFO( logger, "create vector of type " << matrixValueType )

    switch ( matrixValueType )
    {
        case common::scalar::DOUBLE:
            return new DenseVector<double>( distribution, value.getValue<double>() );

        case common::scalar::FLOAT:
            return new DenseVector<float>( distribution, value.getValue<float>() );

        default:
            COMMON_THROWEXCEPTION( "unsupported vector type : " << matrixValueType )
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
    SCAI_ASSERT_DEBUG( localContext, "localContext == NULL" )
    SCAI_ASSERT_DEBUG( haloContext, "haloContext == NULL" )

    // default implementation for matrices that do not support halo context

    if ( localContext.get() != haloContext.get() )
    {
        SCAI_LOG_WARN( logger, *this << ": halo context = " << *haloContext << " ignored" )
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

Matrix& Matrix::operator=( const Matrix& other )
{
    // assignment operator is just implemented by the assign method

    SCAI_LOG_INFO( logger, *this << ": operator = " << other )

    this->assign( other );

    SCAI_LOG_INFO( logger, *this << ": end operator = " << other )

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator=( const Expression_SM& exp )
{
    // exp is Expression object that stands for s * A

    const Matrix& A = exp.getArg2();
    const Scalar& s = exp.getArg1();
    this->matrixTimesScalar( A, s );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator*=( const Scalar exp )
{
    // this *= alpha  -> this->scale( exp )

    this->scale( exp );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator+=( const Expression_SM& exp )
{
    // this += alpha * A  -> this = alpha * A + 1.0 * this

    *this = Expression_SM_SM( exp, Expression_SM( Scalar( 1 ), *this ) );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator-=( const Expression_SM& exp )
{
    // this -= alpha * A  -> this = 1.0 * this + ( - alpha ) * A

    Expression_SM minusExp( -exp.getArg1(), exp.getArg2() );

    *this = Expression_SM_SM( Expression_SM( Scalar( 1 ), *this ), minusExp );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator+=( const Matrix& exp )
{
    // this += A  -> this = 1.0 * A + 1.0 * this

    *this = Expression_SM_SM( Expression_SM( Scalar( 1 ), *this ), Expression_SM( Scalar( 1 ), exp ) );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator-=( const Matrix& exp )
{
    // this -= A  -> this = -1.0 * A + 1.0 * this

    *this = Expression_SM_SM( Expression_SM( Scalar( 1 ), *this ), Expression_SM( Scalar( -1 ), exp ) );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator=( const Expression_SMM& exp )
{
    // exp is Expression object that stands for A * B with matrices A * B
    //   ->   1.0 * A * B + 0.0 * A

    Expression_SM exp2( Scalar( 0 ), *this );

    *this = Expression_SMM_SM( exp, exp2 );

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
    return ( double ) getNumValues() / mNumRows / mNumColumns;
}

/* ---------------------------------------------------------------------------------*/

bool Matrix::checkSymmetry() const
{
    // check symmetry of matrix
    IndexType n = getNumRows();

    if ( n != getNumColumns() )
    {
        return false;
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < i; ++j )
        {
            if ( getValue( i, j ) != getValue( j, i ) )
            {
                return false;
            }
        }
    }

    return true;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::sanityCheck( const Expression<Matrix, Matrix, Times>& exp )
{
    // check sanity of matrix product exp = A * B

    const Matrix& A = exp.getArg1();
    const Matrix& B = exp.getArg2();

    const Distribution& colDistA = A.getColDistribution();
    const Distribution& rowDistB = B.getDistribution();

    if ( colDistA != rowDistB )
    {
        COMMON_THROWEXCEPTION(
            "A * B with A = " << A << ", B = " << B << std::endl << "col size/distribution of A  = " << A.getColDistribution() << " does not match row/size distribution of B = " << B.getDistribution() );
    }
}

void Matrix::sanityCheck( const Expression<Matrix, Matrix, Times>& exp, const Matrix& C )
{
    sanityCheck( exp ); // verify the sanity of the matrix product

    // verify that result of matrix multiplication and C are conform

    const Matrix& A = exp.getArg1();
    const Matrix& B = exp.getArg2();

    const Distribution& rowDistA = A.getDistribution();
    const Distribution& colDistB = B.getColDistribution();

    const Distribution& rowDistC = C.getDistribution();
    const Distribution& colDistC = C.getColDistribution();

    if ( rowDistA != rowDistC )
    {
        COMMON_THROWEXCEPTION( "Size/distribution of rows do not match: " << "ARG1 = " << A << ", ARG2 = " << C )
    }

    if ( colDistB != colDistC )
    {
        COMMON_THROWEXCEPTION( "Size/distribution of cols do not match: " << "ARG1 = " << B << ", ARG2 = " << C )
    }
}

void Matrix::sanityCheck( const Matrix& A, const Matrix& B )
{
    // verify that A and B are conform for addition

    const Distribution& rowDistA = A.getDistribution();
    const Distribution& colDistA = A.getColDistribution();

    const Distribution& rowDistB = B.getDistribution();
    const Distribution& colDistB = B.getColDistribution();

    if ( rowDistA != rowDistB )
    {
        COMMON_THROWEXCEPTION( "Size/distribution of rows do not match: " << "ARG1 = " << A << ", ARG2 = " << B )
    }

    if ( colDistA != colDistB )
    {
        COMMON_THROWEXCEPTION( "Size/distribution of cols do not match: " << "ARG1 = " << A << ", ARG2 = " << B )
    }
}

/* ---------------------------------------------------------------------------------*/

/**
 * @brief the assignment operator for a GEMM expression.
 */
Matrix& Matrix::operator=( const Expression_SMM_SM& exp )
{
    const Expression_SMM& arg1 = exp.getArg1();
    const Expression_SM& arg11 = arg1.getArg1();
    const Expression_SM& arg2 = exp.getArg2();

    const Matrix& A = arg11.getArg2();
    const Matrix& B = arg1.getArg2();
    const Matrix& C = arg2.getArg2();
    const Scalar& alpha = arg11.getArg1();
    const Scalar& beta = arg2.getArg1();

    SCAI_LOG_INFO( logger,
                   "operator=:  " << alpha << " * A * B  + " << beta << " * C" " with A = " << A << ", B = " << B << ", C = " << C )

    const Scalar zero( 0 );

    if ( beta == zero )
    {
        sanityCheck( Expression<Matrix, Matrix, Times>( A, B ) );
    }
    else
    {
        sanityCheck( Expression<Matrix, Matrix, Times>( A, B ), C );
    }

    SCAI_LOG_INFO( logger, "Context of this before matrixTimesMatrix = " << this->getContext() )

    A.matrixTimesMatrix( *this, alpha, B, beta, C );

    SCAI_LOG_INFO( logger, "end operator=:  A * B * alpha + C * beta " )

    SCAI_LOG_INFO( logger, "Context of this after matrixTimesMatrix = " << this->getContext() )

    return *this;
}

/* ---------------------------------------------------------------------------------*/

/**
 * @brief the assignment operator for a MM addition.
 */
Matrix& Matrix::operator=( const Expression_SM_SM& exp )
{
    SCAI_LOG_INFO( logger, "operator=:  A * alpha + B * beta " )

    const Matrix& A = exp.getArg1().getArg2();
    const Matrix& B = exp.getArg2().getArg2();
    const Scalar& alpha = exp.getArg1().getArg1();
    const Scalar& beta = exp.getArg2().getArg1();
    const Scalar zero( 0.0 );

    if ( beta == zero )
    {
        // second summand not needed
        this->matrixTimesScalar( A, alpha );
        return *this;
    }

    if ( alpha == zero )
    {
        // first summand not needed
        this->matrixTimesScalar( B, beta );
        return *this;
    }

    // Do sanity checks

    sanityCheck( A, B );

    this->matrixPlusMatrix( alpha, A, beta, B );

    return *this;
}

} /* end namespace lama */

} /* end namespace scai */
