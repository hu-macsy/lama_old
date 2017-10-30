/**
 * @file _Matrix.cpp
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
 * @brief Implementation of methods for all kind of matrices.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/matrix/_Matrix.hpp>

// local library
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/io/PartitionIO.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

namespace scai
{

using namespace common;
using namespace dmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( Matrix::logger, "Matrix" )

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a matrix                                                          */
/* ---------------------------------------------------------------------------------------*/

Matrix* Matrix::getMatrix( const Format::MatrixStorageFormat format, const common::scalar::ScalarType valueType )
{
    MatrixCreateKeyType mattype( format, valueType );
    return Matrix::create( mattype );
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( const Matrix& other ) :

    Distributed( other ),
    mColDistribution( other.mColDistribution ),
    mCommunicationKind( other.mCommunicationKind )
{
    SCAI_LOG_INFO( logger, "Creating copy of " << other << " with same distributions." )
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist ) :

    Distributed( rowDist ),
    mColDistribution( colDist ),
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
    mColDistribution( DistributionPtr( new NoDistribution( numColumns ) ) )
{
    setDefaultKind();
    SCAI_LOG_INFO( logger, "Creating a replicated Matrix of size " << numRows << " x " << numColumns )
}

/* ----------------------------------------------------------------------- */

void Matrix::setIdentity( const IndexType n )
{
    // take replicated distribution and use pure method
    setIdentity( DistributionPtr( new NoDistribution( n ) ) );
}

/* ----------------------------------------------------------------------- */

void Matrix::setDiagonalProperty()
{
    SCAI_ASSERT_EQ_ERROR( getRowDistribution(), getColDistribution(),
                          "col/row distribution must be equal to set diagonal property" );

    // Now we can set it for the local storage

    _MatrixStorage& m = const_cast<_MatrixStorage&>( getLocalStorage() );

    bool errorFlag = false;

    try
    {
        m.setDiagonalProperty();
    }
    catch ( Exception& e )
    {
        SCAI_LOG_ERROR( logger, "This processor could not force diagonal property" )
        errorFlag = true;
    }

    errorFlag = getRowDistribution().getCommunicator().any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Not all processes could set diagonal property" )
    }
}

/* ----------------------------------------------------------------------- */

void Matrix::checkSettings() const
{
    if ( !mColDistribution )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for column distribution" )
    }
}

/* ----------------------------------------------------------------------- */

Matrix::Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
    : Distributed( rowDistribution )
{
    setDistributedMatrix( rowDistribution, colDistribution );
    setDefaultKind();
    SCAI_LOG_INFO( logger,
                   "Construct a Matrix of size " << getNumRows() << " x " << getNumColumns() 
                    << " with the distribution " << getDistribution() )
}

Matrix::Matrix( DistributionPtr distribution )
    : Distributed( distribution )
{
    setDistributedMatrix( distribution, distribution );
    SCAI_LOG_INFO( logger,
                   "Construct a square Matrix of size " << getNumRows() << " x " << getNumColumns() 
                   << " with the row/col distribution " << getDistribution() )
}

Matrix::Matrix() : 

    Distributed( DistributionPtr( new NoDistribution( 0 ) ) ), 
    mColDistribution( DistributionPtr( new NoDistribution( 0 ) ) )
{
    setDefaultKind();
}

Matrix::~Matrix()
{
    SCAI_LOG_DEBUG( logger, "~Matrix" )
}

Matrix::SyncKind Matrix::getDefaultSyncKind()
{
    static bool computed = false;

    static SyncKind syncKind = ASYNCHRONOUS;

    if ( !computed )
    {
        bool isAsync = true;

        common::Settings::getEnvironment( isAsync, "SCAI_ASYNCHRONOUS" );

        if ( !isAsync )
        {
            syncKind = SYNCHRONOUS;
        }
    }

    return syncKind;
}

void Matrix::setDefaultKind()
{
    mCommunicationKind = getDefaultSyncKind();
}

/* ---------------------------------------------------------------------------------*/

void Matrix::buildCSRGraph( IndexType ia[], IndexType ja[], IndexType vwgt[], const IndexType* globalIndexes ) const
{
    getLocalStorage().buildCSRGraph( ia, ja, vwgt, globalIndexes );
}

/* ---------------------------------------------------------------------------------*/

IndexType Matrix::getCSRGraphSize() const
{
    // Currently only supported if column distribution is replicated
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), getLocalStorage().getNumColumns(), "getCSRGraphSize only for replicated column distribution" )
    // diagonal elements will not be used
    return getLocalNumValues() - getDistribution().getLocalSize();
}

/* ---------------------------------------------------------------------------------*/

void Matrix::setDistributedMatrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_ASSERT_ERROR( rowDistribution, "NULL row distribution for matrix not allowed" )
    SCAI_ASSERT_ERROR( colDistribution, "NULL column distribution for matrix not allowed" )
    setDistributionPtr( rowDistribution );
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

void Matrix::inheritAttributes( const Matrix& other )
{
    setCommunicationKind( other.getCommunicationKind() );
    setContextPtr( other.getContextPtr() );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::writeAt( std::ostream& stream ) const
{
    stream << "Matrix(" << getNumRows() << "x" << getNumColumns() << ")";
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
    *this = Expression_SM_SM( exp, Expression_SM( Scalar( 1.0 ), *this ) );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator-=( const Expression_SM& exp )
{
    // this -= alpha * A  -> this = 1.0 * this + ( - alpha ) * A
    Expression_SM minusExp( -exp.getArg1(), exp.getArg2() );
    *this = Expression_SM_SM( Expression_SM( Scalar( 1.0 ), *this ), minusExp );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator+=( const Matrix& exp )
{
    // this += A  -> this = 1.0 * A + 1.0 * this
    *this = Expression_SM_SM( Expression_SM( Scalar( 1.0 ), *this ), Expression_SM( Scalar( 1.0 ), exp ) );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator-=( const Matrix& exp )
{
    // this -= A  -> this = -1.0 * A + 1.0 * this
    *this = Expression_SM_SM( Expression_SM( Scalar( 1.0 ), *this ), Expression_SM( Scalar( -1.0 ), exp ) );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

Matrix& Matrix::operator=( const Expression_SMM& exp )
{
    // exp is Expression object that stands for A * B with matrices A * B
    //   ->   1.0 * A * B + 0.0 * A
    Expression_SM exp2( Scalar( 0.0 ), *this );
    *this = Expression_SMM_SM( exp, exp2 );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::swapMatrix( Matrix& other )
{
    Distributed::swap( other );
    std::swap( mColDistribution, other.mColDistribution );
}

/* ---------------------------------------------------------------------------------*/

double Matrix::getSparsityRate() const
{
    return ( double ) getNumValues() / getNumRows() / getNumColumns();
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

    // Note: this solution is not very efficient

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

    SCAI_LOG_INFO( logger, "Context of this before matrixTimesMatrix = " << *getContextPtr() )
    A.matrixTimesMatrix( *this, alpha, B, beta, C );
    SCAI_LOG_INFO( logger, "end operator=:  A * B * alpha + C * beta " )
    SCAI_LOG_INFO( logger, "Context of this after matrixTimesMatrix = " << *getContextPtr() )
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
        // second term not needed
        this->matrixTimesScalar( A, alpha );
        return *this;
    }

    if ( alpha == zero )
    {
        // first term not needed
        this->matrixTimesScalar( B, beta );
        return *this;
    }

    // Do sanity checks
    sanityCheck( A, B );
    this->matrixPlusMatrix( alpha, A, beta, B );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::writeToSingleFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType /* = UNKNOWN for DEFAULT */,
    const common::scalar::ScalarType indexType /* = UNKNOWN for DEFAULT */,
    const FileIO::FileMode fileMode /* = DEFAULT_MODE */ ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    if ( getDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        if ( comm->getRank() == 0 )
        {
            getLocalStorage().writeToFile( fileName, fileType, dataType, indexType, fileMode );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file
        comm->synchronize();
    }
    else
    {
        DistributionPtr rowDist( new NoDistribution( getNumRows() ) );
        DistributionPtr colDist( new NoDistribution( getNumColumns() ) );
        common::unique_ptr<Matrix> repM( copy( rowDist, colDist ) );
        repM->writeToSingleFile( fileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::writeToPartitionedFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType /* = UNKNOWN for DEFAULT */,
    const common::scalar::ScalarType indexType /* = UNKNOWN for DEFAULT */,
    const FileIO::FileMode fileMode /* = DEFAULT_MODE */ ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    if ( getColDistribution().isReplicated() )
    {
        // each processor writes its partition to a file with unique name

        getLocalStorage().writeToFile( fileName, fileType, dataType, indexType, fileMode );
    }
    else
    {
        DistributionPtr colDist( new NoDistribution( getNumColumns() ) );
        common::unique_ptr<Matrix> repM( copy( getRowDistributionPtr(), colDist ) );
        repM->writeToPartitionedFile( fileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::writeToFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType,
    const common::scalar::ScalarType indexType,
    const FileIO::FileMode fileMode ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    std::string newFileName = fileName;

    bool isPartitioned;

    const Communicator& comm = getRowDistribution().getCommunicator();

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, comm );

    if ( !isPartitioned )
    {
        writeToSingleFile( newFileName, fileType, dataType, indexType, fileMode );
    }
    else
    {
        writeToPartitionedFile( newFileName, fileType, dataType, indexType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::readFromSingleFile( const std::string& fileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const PartitionId MASTER = 0;
    const PartitionId myRank = comm->getRank();

    // this is a bit tricky stuff, but it avoids an additional copy from storage -> matrix

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    IndexType dims[2];

    if ( myRank == MASTER )
    {
        localMatrix.readFromFile( fileName );

        dims[0] = localMatrix.getNumRows();
        dims[1] = localMatrix.getNumColumns();
    }

    comm->bcast( dims, 2, MASTER );

    if ( myRank != MASTER )
    {
        IndexType localNumRows = 0;
        localMatrix.allocate( localNumRows, dims[1] );
    }

    DistributionPtr rowDist( new SingleDistribution( dims[0], comm, MASTER ) );
    DistributionPtr colDist( new NoDistribution( dims[1] ) );

    // works fine as assign can deal with alias, i.e. localMatrix und getLocalStorage() are same

    SCAI_LOG_DEBUG( logger, *comm << ": assign local storage " << localMatrix );

    assign( localMatrix, rowDist, colDist );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::readFromSingleFile( const std::string& fileName, const DistributionPtr distribution )
{
    if ( distribution.get() == NULL )
    {
        readFromSingleFile( fileName );
        return;
    }

    // dist must be block distributed, not checked again here

    const IndexType n = distribution->getBlockDistributionSize();

    if ( n == nIndex )
    {
        readFromSingleFile( fileName );
        redistribute( distribution, getColDistributionPtr() );
        return;
    }

    const Communicator& comm = distribution->getCommunicator();

    IndexType first = 0;

    if ( n > 0 )
    {
        first = distribution->local2global( 0 );   // first global index
    }

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    bool error = false;

    try
    {
        localMatrix.readFromFile( fileName, first, n );
    }
    catch ( Exception& ex )
    {
        SCAI_LOG_ERROR( logger, ex.what() )
        error = true;
    }

    error = distribution->getCommunicator().any( error );

    if ( error )
    {
        COMMON_THROWEXCEPTION( "readFromSingleFile failed." )
    }

    IndexType numColumns = comm.max( localMatrix.getNumColumns() );

    DistributionPtr colDist( new NoDistribution( numColumns ) );

    assign( localMatrix, distribution, colDist );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::readFromPartitionedFile( const std::string& myPartitionFileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // this is a bit tricky stuff, but it avoids an additional copy from storage -> matrix

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( getLocalStorage() );

    bool errorFlag = false;

    IndexType localSize = 0;

    try
    {
        localMatrix.readFromFile( myPartitionFileName );

        localSize = localMatrix.getNumRows();

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

    // We assume a general block distribution

    IndexType globalSize = comm->sum( localSize );

    DistributionPtr rowDist( new GenBlockDistribution( globalSize, localSize, comm ) );

    // make sure that all processors have the same number of columns

    IndexType numColumns = comm->max( localMatrix.getNumColumns() );

    // for consistency we have to set the number of columns in each stroage

    localMatrix.setDimension( localSize, numColumns );

    DistributionPtr colDist( new NoDistribution( numColumns ) );

    assign( localMatrix, rowDist, colDist );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::resetRowDistribution( DistributionPtr newDist )
{
    SCAI_ASSERT_EQ_ERROR( getNumRows(), newDist->getGlobalSize(), "global size mismatch" )

    const _MatrixStorage& localMatrix = getLocalStorage();

    SCAI_ASSERT_EQ_ERROR( localMatrix.getNumRows(), newDist->getLocalSize(), "local size mismatch" );

    setDistributionPtr( newDist );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::resetRowDistributionByFirstColumn()
{
    if ( getRowDistribution().isReplicated() )
    {
        return;   // nothing to do
    }

    bool errorFlag = false;

    CommunicatorPtr comm = getRowDistribution().getCommunicatorPtr();

    // catch local exceptions and throw later a global exception

    try
    {
        SCAI_LOG_INFO( logger, "getRowDistributionByFirstColumn" )

        const _MatrixStorage& localMatrix = getLocalStorage();

        hmemo::HArray<IndexType> myGlobalIndexes;

        localMatrix.getFirstColumnIndexes( myGlobalIndexes );

        SCAI_LOG_DEBUG( logger, "first col indexes = " << myGlobalIndexes )

        // if storage has not the global column index of diagonal first, this test is likely to fail

        SCAI_ASSERT_DEBUG( utilskernel::HArrayUtils::isSorted( myGlobalIndexes, common::binary::LE ),
                           "first column indexes are not sorted, cannot be global indexes" )

        // otherwise building the distribution will fail

        DistributionPtr dist( new dmemo::GeneralDistribution( getNumRows(), myGlobalIndexes, comm ) );

        resetRowDistribution( dist );
    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_ERROR( logger, *comm << ": serious error for building general distribution by first col index" )
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "determing general distribution by column indexes failed." )
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::readFromFile( const std::string& matrixFileName, const std::string& distributionFileName )
{
    if ( distributionFileName.size() == 0 )
    {
        readFromFile( matrixFileName );
        resetRowDistributionByFirstColumn();
    }
    else if ( distributionFileName == "BLOCK" )
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        DistributionPtr rowDist;

        // for a single file we set a BlockDistribution

        if ( matrixFileName.find( "%r" ) == std::string::npos )
        {
            PartitionId root = 0;

            IndexType numRows = nIndex;

            if ( comm->getRank() == root )
            {
                numRows = FileIO::getStorageSize( matrixFileName );
            }

            comm->bcast( &numRows, 1, root );

            rowDist.reset( new BlockDistribution( numRows, comm ) );
        }

        readFromFile( matrixFileName, rowDist );
    }
    else
    {
        // read the distribution

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        DistributionPtr rowDist = PartitionIO::readDistribution( distributionFileName, comm );

        readFromFile( matrixFileName, rowDist );
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::readFromFile( const std::string& fileName, DistributionPtr rowDist )
{
    SCAI_LOG_INFO( logger,
                   *this << ": readFromFile( " << fileName << " )" )

    std::string newFileName = fileName;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // take default

    if ( rowDist.get() )
    {
        comm = rowDist->getCommunicatorPtr();
    }

    bool isPartitioned;

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, *comm );

    SCAI_LOG_INFO( logger, *comm << ": Matrix.readFromFile ( " << fileName << " ) -> read "
                   << newFileName << ", partitioned = " << isPartitioned );

    if ( !isPartitioned )
    {
        readFromSingleFile( newFileName, rowDist );
    }
    else
    {
        readFromPartitionedFile( newFileName );

        if ( rowDist.get() )
        {
            resetRowDistribution( rowDist );
        }
    }
}

/* ---------------------------------------------------------------------------------*/

void Matrix::concatenate( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist, const std::vector<const Matrix*>& matrices )
{
    COMMON_THROWEXCEPTION( "concatenation of matrices not supported, #matrices = " << matrices.size()
                           << ", row dist = " << *rowDist << ", col dist = " << *colDist )
}

/* ---------------------------------------------------------------------------------*/

void Matrix::vcat( const Matrix& m1, const Matrix& m2 )
{
    SCAI_ASSERT_EQ_ERROR( m1.getRowDistribution(), m2.getRowDistribution(), "vcat: matrices must have same row distribution" )

    DistributionPtr rowDist = m1.getRowDistributionPtr();

    DistributionPtr colDist( new NoDistribution( m1.getNumColumns() + m2.getNumColumns() ) );

    std::vector<const Matrix*> matrices;

    matrices.push_back( &m1 );
    matrices.push_back( &m2 );
    
    concatenate( rowDist, colDist, matrices );
}

/* ---------------------------------------------------------------------------------*/

void Matrix::hcat( const Matrix& m1, const Matrix& m2 )
{
    SCAI_ASSERT_EQ_ERROR( m1.getNumColumns(), m2.getNumColumns(), "No horizontal cut possible due to different column sizes" )
 
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr rowDist( new BlockDistribution( m1.getNumRows() + m2.getNumRows(), comm ) ); 
    DistributionPtr colDist( new NoDistribution( m1.getNumColumns() ) );
    
    std::vector<const Matrix*> matrices;

    matrices.push_back( &m1 );
    matrices.push_back( &m2 );

    concatenate( rowDist, colDist, matrices );
}

/* ---------------------------------------------------------------------------------*/

Scalar Matrix::maxDiffNorm( const Matrix& other ) const
{
    IndexType nRows = getNumRows();
    IndexType nCols = getNumColumns();
    SCAI_ASSERT_EQUAL( nRows, other.getNumRows(), "size mismatch" )
    SCAI_ASSERT_EQUAL( nCols, other.getNumColumns(), "size mismatch" )
    VectorCreateKeyType vectorType1( _Vector::DENSE, getValueType() );
    VectorCreateKeyType vectorType2( _Vector::DENSE, other.getValueType() );
    common::unique_ptr<_Vector> ptrRow1( _Vector::create( vectorType1 ) );
    common::unique_ptr<_Vector> ptrRow2( _Vector::create( vectorType2 ) );
    Scalar diff( 0 );

    // now traverse  all rows

    for ( IndexType i = 0; i < nRows; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices
        getRow( *ptrRow1, i );
        other.getRow( *ptrRow2, i );

        // compare the two vectors element-wise

        for ( IndexType j = 0; j < nCols; j++ )
        {
            Scalar elem1 = ptrRow1->getValue( j );
            Scalar elem2 = ptrRow2->getValue( j );
            Scalar diff1  = abs( elem1 - elem2 );

            if ( diff1 > diff )
            {
                diff = diff1;
            }
        }
    }

    return diff;
}

/* ---------------------------------------------------------------------------------*/

void Matrix::redistribute( const dmemo::Redistributor& redistributor )
{
    if ( getColDistribution().isReplicated() ) 
    {
        redistribute( redistributor, getColDistributionPtr() );
    }
    else if ( getColDistribution() == getRowDistribution() )
    {
        redistribute( redistributor, redistributor.getTargetDistributionPtr() );
    }
    else
    {
        COMMON_THROWEXCEPTION( "redistribute: no new column distribution" )
    }
}

/* ---------------------------------------------------------------------------------*/

Matrix* Matrix::copy( DistributionPtr rowDistribution, DistributionPtr colDistribution ) const
{
    // simple default implementation that works for each matrix
    common::unique_ptr<Matrix> rep( copy() );
    // unique_ptr guarantees that data is freed if redistribute fails for any reason
    rep->redistribute( rowDistribution, colDistribution );
    return rep.release();
}

MatrixCreateKeyType Matrix::getCreateValue() const
{
    return MatrixCreateKeyType( getFormat(), getValueType() );
}

} /* end namespace lama */

} /* end namespace scai */
