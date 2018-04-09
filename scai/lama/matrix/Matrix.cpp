/**
 * @file Matrix.cpp
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
 * @brief Implementation of methods for the abstract class Matrix<ValueType>
 * @author Thomas Brandes
 * @date 31.10.2017
 */

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/tracing.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/common/macros/unsupported.hpp>
#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using common::TypeTraits;
using namespace dmemo;
using hmemo::HArray;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*    static methods                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Matrix<ValueType>* Matrix<ValueType>::getMatrix( Format format )
{
    return static_cast<Matrix<ValueType>*>( _Matrix::getMatrix( format, TypeTraits<ValueType>::stype ) );
}

/* ------------------------------------------------------------------------- */
/*    Constructors / Destructor                                              */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Matrix<ValueType>::Matrix() : _Matrix()
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const IndexType numRows, const IndexType numColumns ) : 

    _Matrix( numRows, numColumns )

{
    SCAI_LOG_DEBUG( logger, "Matrix<" << TypeTraits<ValueType>::id() << "> ( "
                            << _Matrix::getNumRows() << " x " << _Matrix::getNumColumns() << " )" )
}

template<typename ValueType>
Matrix<ValueType>::Matrix( DistributionPtr rowDistribution, DistributionPtr colDistribution ) :

    _Matrix( rowDistribution, colDistribution )
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const _Matrix& other, 
                           DistributionPtr rowDistribution, 
                           DistributionPtr colDistribution ) :

    _Matrix( other, rowDistribution, colDistribution )

{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const _Matrix& other ) :

    _Matrix( other )
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const Matrix<ValueType>& other ) :

    _Matrix( other )
{
}

template<typename ValueType>
Matrix<ValueType>::~Matrix()
{
    SCAI_LOG_DEBUG( logger, "~Matrix<" << TypeTraits<ValueType>::id() << ">" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType Matrix<ValueType>::operator()( IndexType i, IndexType j ) const
{
    return getValue( i, j );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType Matrix<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
size_t Matrix<ValueType>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

/* ========================================================================= */

template<typename ValueType>
bool Matrix<ValueType>::checkSymmetry() const
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

/* ========================================================================= */

template<typename ValueType>
Vector<ValueType>* Matrix<ValueType>::newTargetVector() const
{
    std::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>( getContextPtr() ) );
    vector->setSameValue( this->getRowDistributionPtr(), 0 );
    return vector.release();
}

/* ========================================================================= */

template<typename ValueType>
Vector<ValueType>* Matrix<ValueType>::newSourceVector() const
{
    std::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>( getContextPtr() ) );
    vector->setSameValue( this->getColDistributionPtr(), 0 );
    return vector.release();
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::matrixTimesVector(
    Vector<ValueType>& result,
    const ValueType alpha,
    const Vector<ValueType>& x,
    const ValueType beta,
    const Vector<ValueType>* y,
    common::MatrixOp op ) const
{
    SCAI_REGION( "Mat.timesVector" )

    SCAI_LOG_INFO( logger, 
                   "result = " << alpha << " * M<" << this->getValueType() << ","
                   << op 
                   << ">[" << this->getNumRows() << " x " << this->getNumColumns() << "]"
                   << " * x[" << x.size() << "] + " << beta << " * y[]" )

    if ( beta == common::Constants::ZERO )
    {
        if ( y != nullptr )
        {
            SCAI_LOG_INFO( logger, "this vector is ignored (beta == 0) : " << y )
            matrixTimesVector( result, alpha, x, beta, nullptr, op );
            return;
        }
    }
    else
    {
        SCAI_ASSERT_ERROR( y != nullptr, "vector y is null pointer, but beta != 0" )
    }

    DistributionPtr sourceDist = common::isTranspose( op ) ? getRowDistributionPtr() : getColDistributionPtr();
    DistributionPtr targetDist = common::isTranspose( op ) ? getColDistributionPtr() : getRowDistributionPtr();

    // temorary X required if not DENSE, distribution does not match or if an alias

    bool needsTemporaryX = false;

    if ( &x == &result )
    {
        SCAI_UNSUPPORTED( "z = alpha * A * x + beta * y: temporary needed for x due to alias with z" );
        needsTemporaryX = true;
    }

    if ( &x == y ) 
    {
        SCAI_UNSUPPORTED( "z = alpha * A * x + beta * y: temporary needed for x due to alias with y" );
        needsTemporaryX = true;
    }

    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "z = alpha * A * x + beta * y: temporary needed for x as it is sparse" );
        needsTemporaryX = true;
    }

    if ( x.getDistribution() != *sourceDist )
    {
        SCAI_UNSUPPORTED( "z = alpha * A * x + beta * y: temporary needed for x as distribution does not match" );
        needsTemporaryX = true;
    }

    if ( needsTemporaryX )
    {
        matrixTimesVector( result, alpha, distribute<DenseVector<ValueType>>( x, sourceDist ), beta, y, op );  
        return;
    }

    bool needsTemporaryY = false;

    if ( y != nullptr )
    {
        if ( y->getVectorKind() != VectorKind::DENSE )
        {
            SCAI_UNSUPPORTED( "matrixTimesVector: temporary needed for y as it is sparse vector" );
            needsTemporaryY = true;
        }
    
        if ( y->getDistribution() != *targetDist )
        {
            SCAI_UNSUPPORTED( "matrixTimesVector: temporary needed for y as distribution does not match" );
            needsTemporaryY = true;
        }
    }

    if ( needsTemporaryY )
    {
        auto tmpY = distribute<DenseVector<ValueType>>( *y, targetDist );
        matrixTimesVector( result, alpha, x, beta, &tmpY, op );
        return;
    }

    if ( result.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "matrixTimesVector: temporary needed for result as not dense" )
        DenseVector<ValueType> tmpResult;
        tmpResult.allocate( targetDist );            // no initialization required
        matrixTimesVector( tmpResult, alpha, x, beta, y, op );
        result = tmpResult;
        return;
    }

    if ( &result == y  )
    {
        SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
    }
    else
    {
        // we inherit the row distribution of this matrix to result vector

        result.allocate( targetDist );
    }

    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>* denseY = nullptr;   // optional

    if ( y != nullptr )
    {
        denseY = static_cast<const DenseVector<ValueType>*>( y );
    }

    DenseVector<ValueType>& denseResult = static_cast<DenseVector<ValueType>&>( result );

    // Now call the version with dense vector implemented by derived class

    matrixTimesVectorDense( denseResult, alpha, denseX, beta, denseY, op );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::setRow( 
    const Vector<ValueType>& row, 
    const IndexType globalRowIndex,
    const common::BinaryOp op )
{
    SCAI_ASSERT_EQ_ERROR( row.size(), this->getNumColumns(), "row size mismatch" )

    SCAI_LOG_DEBUG( logger, "setRow " << globalRowIndex << ": row = " << row << ", op = " << op )

    bool needsTmp = false;

    if ( ! row.getDistribution().isReplicated() )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setRow, row is not replicated, use temporary" )
    }
    if ( row.getVectorKind() != VectorKind::DENSE )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setRow, row is not DENSE vector" )
    }

    if ( needsTmp )
    {
        DenseVector<ValueType> tmpRow;
        tmpRow.assign( row );
        tmpRow.replicate();
        setRow( tmpRow, globalRowIndex, op );
        return;
    }

    SCAI_ASSERT_VALID_INDEX_ERROR( globalRowIndex, this->getNumRows(), "illegal row index" )

    // row should be a DenseVector of same type, otherwise use a temporary

    const DenseVector<ValueType>& denseRow = static_cast<const DenseVector<ValueType>&>( row );

    // owner sets the row, maybe each processor for replicated row distribution

    IndexType localRowIndex = this->getRowDistribution().global2local( globalRowIndex );

    if ( localRowIndex != invalidIndex )
    {
        this->setLocalRow( denseRow.getLocalValues(), localRowIndex, op );
    }
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::setColumn( 
    const Vector<ValueType>& column,
    const IndexType colIndex,
    const common::BinaryOp op )
{
    using namespace scai::hmemo;

    bool needsTmp = false;

    if ( column.getDistribution() != this->getRowDistribution() )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setColumn, distribution of col does not match distribution of matrix, use temporary" )
    }
    if ( column.getVectorKind() != VectorKind::DENSE )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setColumn, col is not DENSE vector" )
    }

    if ( needsTmp )
    {
        setColumn( distribute<DenseVector<ValueType>>( column, getRowDistributionPtr() ), colIndex, op );
        return;
    }

    SCAI_ASSERT_VALID_INDEX_ERROR( colIndex, this->getNumColumns(), "illegal col index" )

    const DenseVector<ValueType>& denseColumn = static_cast<const DenseVector<ValueType>&>( column );

    SCAI_ASSERT_EQ_ERROR( denseColumn.getDistribution(), this->getRowDistribution(), 
                          "distribution of vector to set as column must match distribution of rows" )

    this->setLocalColumn( denseColumn.getLocalValues(), colIndex, op );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::vectorTimesMatrixRepCols(
    DenseVector<ValueType>& denseResult,
    const ValueType alpha,
    const DenseVector<ValueType>& denseX,
    const ValueType beta,
    const DenseVector<ValueType>* denseY ) const
{
    SCAI_REGION( "Mat.vectorTimesMatrixRepCols" )

    hmemo::HArray<ValueType>& localResult = denseResult.getLocalValues();

    // be careful: denseY is undefined if beta == 0

    const HArray<ValueType>& localY = denseY == nullptr ? localResult : denseY->getLocalValues();
    const HArray<ValueType>& localX = denseX.getLocalValues();

    const Distribution& colDist = this->getColDistribution();

    // this routine is only for non-replicated columns, i.e. mHaloData is empty

    SCAI_ASSERT( 1, colDist.getNumPartitions() );

    const Distribution& rowDist = this->getRowDistribution();
    const Communicator& comm = rowDist.getCommunicator();

    const MatrixStorage<ValueType>& localData = this->getLocalStorage();

    if ( comm.getRank() == 0 )
    {
        // only one single processor adds beta * y
        localData.matrixTimesVector( localResult, alpha, localX, beta, localY, common::MatrixOp::TRANSPOSE );
    }
    else
    {
        localData.matrixTimesVector( localResult, alpha, localX, ValueType( 0 ), localY, common::MatrixOp::TRANSPOSE );
    }

    if ( comm.getSize() >  1 )
    {
        // Sum up all incarnations of localResult

        comm.sumArray( localResult );
    }
}

/* ========================================================================= */

template<typename ValueType>
RealType<ValueType> Matrix<ValueType>::maxDiffNorm( const Matrix<ValueType>& other ) const
{
    IndexType nRows = getNumRows();
    IndexType nCols = getNumColumns();

    SCAI_ASSERT_EQUAL( nRows, other.getNumRows(), "size mismatch" )
    SCAI_ASSERT_EQUAL( nCols, other.getNumColumns(), "size mismatch" )

    DenseVector<ValueType> row;
    DenseVector<ValueType> rowOther;

    RealType<ValueType> diff = 0;

    // now traverse  all rows

    for ( IndexType i = 0; i < nRows; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices

        getRow( row, i );
        other.getRow( rowOther, i );

        RealType<ValueType> diffRow = row.maxDiffNorm( rowOther );

        if ( diffRow > diff )
        {
            diff = diffRow;
        }
    }

    return diff;
}

/* ========================================================================= */

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const Matrix<ValueType>& other )
{
    this->assign( other );
    return *this;
}

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const OpMatrix<ValueType>& other )
{
    const Matrix<ValueType>& matrix = other.getMatrix();
    const common::MatrixOp op = other.getOp();

    switch ( op ) 
    {
        case common::MatrixOp::NORMAL:  
            this->assign( matrix ); 
            break;
        case common::MatrixOp::TRANSPOSE:  
            this->assignTranspose( matrix ); 
            break;
        default:
            COMMON_THROWEXCEPTION( "matrix = " << op << "( matrix ) not supported yet" )
    }

    return *this;
}

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const Expression_SMM_SM<ValueType>& exp )
{
    const Expression_SMM<ValueType>& arg1 = exp.getArg1();
    const Expression_SM<ValueType>& arg11 = arg1.getArg1();
    const Expression_SM<ValueType>& arg2 = exp.getArg2();

    const OpMatrix<ValueType>& opMatA = arg11.getArg2();
    const OpMatrix<ValueType>& opMatB = arg1.getArg2();
    const OpMatrix<ValueType>& opMatC = arg2.getArg2();

    const Matrix<ValueType>& matA = opMatA.getMatrix();
    const Matrix<ValueType>& matB = opMatB.getMatrix();
    const Matrix<ValueType>& matC = opMatC.getMatrix();

    const intern::Scalar& alphaS = arg11.getArg1();
    const intern::Scalar& betaS  = arg2.getArg1();

    const ValueType& alpha = alphaS.getValue<ValueType>();
    const ValueType& beta  = betaS.getValue<ValueType>();

    SCAI_LOG_INFO( logger,
                   "operator=:  " << alpha << " * A * B  + " << beta << " * C"
                   << " with A = " << matA << ", B = " << matB << ", C = " << matC )

    SCAI_ASSERT_EQ_ERROR( opMatA.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )
    SCAI_ASSERT_EQ_ERROR( opMatB.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )
    SCAI_ASSERT_EQ_ERROR( opMatC.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )

    matA.matrixTimesMatrix( *this, alpha, matB, beta, matC );

    SCAI_LOG_INFO( logger, "Context of this after matrixTimesMatrix = " << *getContextPtr() )

    return *this;
}


template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const Expression_SM_SM<ValueType>& exp )
{
    SCAI_LOG_INFO( logger, "operator=:  A * alpha + B * beta " )

    const OpMatrix<ValueType>& opMatA = exp.getArg1().getArg2();
    const OpMatrix<ValueType>& opMatB = exp.getArg2().getArg2();

    const Matrix<ValueType>& matA = opMatA.getMatrix();
    const Matrix<ValueType>& matB = opMatB.getMatrix();

    SCAI_ASSERT_EQ_ERROR( opMatA.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )
    SCAI_ASSERT_EQ_ERROR( opMatB.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )

    const intern::Scalar& alphaS = exp.getArg1().getArg1();
    const intern::Scalar& betaS = exp.getArg2().getArg1();

    const ValueType& alpha = alphaS.getValue<ValueType>();
    const ValueType& beta = betaS.getValue<ValueType>();

    const ValueType zero = 0;

    if ( beta == zero )
    {
        // second term not needed
        this->matrixTimesScalar( matA, alpha );
        return *this;
    }

    if ( alpha == zero )
    {
        // first term not needed
        this->matrixTimesScalar( matB, beta );
        return *this;
    }

    // conformance check of matrices A and B is done in the routines

    this->matrixPlusMatrix( alpha, matA, beta, matB );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const Expression_SMM<ValueType>& exp )
{
    // exp =  ( alpha * A ) * B

    const intern::Scalar& alpha = exp.getArg1().getArg1();

    const OpMatrix<ValueType>& opMatA = exp.getArg1().getArg2();
    const OpMatrix<ValueType>& opMatB = exp.getArg2();

    const Matrix<ValueType>& matA = opMatA.getMatrix();
    const Matrix<ValueType>& matB = opMatB.getMatrix();

    SCAI_ASSERT_EQ_ERROR( opMatA.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )
    SCAI_ASSERT_EQ_ERROR( opMatB.getOp(), common::MatrixOp::NORMAL, "unsupported exp" )

    matA.matrixTimesMatrix( *this, alpha.getValue<ValueType>(), matB, ValueType( 0 ), *this );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator=( const Expression_SM<ValueType>& exp )
{
    // exp is Expression object that stands for s * A

    const OpMatrix<ValueType>& opA = exp.getArg2();

    common::MatrixOp op = opA.getOp();
    SCAI_ASSERT_EQ_ERROR( op, common::MatrixOp::NORMAL, "matrix op = " << op << " unsupported in matrixA = alpha * op( matrixB )" )
    const intern::Scalar& s = exp.getArg1();
    this->matrixTimesScalar( opA.getMatrix(), s.getValue<ValueType>() );
    return *this;
}


/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator+=( const Matrix<ValueType>& mat )
{
    // this += A  -> this = 1.0 * A + 1.0 * this

    this->matrixPlusMatrix( ValueType( 1 ), *this, ValueType( 1 ), mat );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator-=( const Matrix<ValueType>& mat )
{
    // this -= A  -> this = -1.0 * A + 1.0 * this

    this->matrixPlusMatrix( ValueType( 1 ), *this, ValueType( -1 ), mat );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator+=( const Expression_SM<ValueType>& exp )
{
    // this += alpha * A  -> this = alpha * A + 1.0 * this

    const OpMatrix<ValueType>& opMat = exp.getArg2();

    SCAI_ASSERT_EQ_ERROR( opMat.getOp(), common::MatrixOp::NORMAL, "unsupported matrix op" )

    const intern::Scalar& s = exp.getArg1();

    this->matrixPlusMatrix( ValueType( 1 ), *this, s.getValue<ValueType>(), opMat.getMatrix() );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator-=( const Expression_SM<ValueType>& exp )
{
    // this -= alpha * A  -> this = 1.0 * this + ( - alpha ) * A

    const OpMatrix<ValueType>& opMat = exp.getArg2();

    SCAI_ASSERT_EQ_ERROR( opMat.getOp(), common::MatrixOp::NORMAL, "unsupported matrix op" )

    const intern::Scalar& s = exp.getArg1();

    this->matrixPlusMatrix( ValueType( 1 ), *this, -s.getValue<ValueType>(), opMat.getMatrix() );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Matrix<ValueType>& Matrix<ValueType>::operator*=( const ValueType alpha )
{
    this->scale( alpha );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void Matrix<ValueType>::concatenate( 
    DistributionPtr rowDist, 
    DistributionPtr colDist, 
    const std::vector<const Matrix<ValueType>*>& matrices )
{
    SCAI_LOG_INFO( logger, "Concatenate " << matrices.size() << " matrices into  a new matrix." )

    // Each processor assembles the local part of each input matrix for the result matrix

    MatrixAssembly<ValueType> assembly;

    IndexType offsetRow = 0;    // row offset where input matrix is set in result matrix
    IndexType offsetCol = 0;    // col offset where input matrix is set in result matrix

    for ( size_t k = 0; k < matrices.size(); ++k )
    {
        const Matrix<ValueType>& m = *matrices[k];

        SCAI_LOG_DEBUG( logger, "dissassemble this matrix: " << m )

        if ( offsetRow + m.getNumRows() > rowDist->getGlobalSize() )
        {
            COMMON_THROWEXCEPTION( "concatenation fails, this arg fails: " << m )
        }

        if ( offsetCol + m.getNumColumns() > colDist->getGlobalSize() )
        {
            COMMON_THROWEXCEPTION( "concatenation fails, arg " << k << " fails: " << m )
        }

        m.disassemble( assembly, offsetRow, offsetCol );

        offsetCol += m.getNumColumns();

        // decide by the sizes where (horizontally or vertically)  we concatenate the next 

        if ( offsetCol == colDist->getGlobalSize() )
        {
            offsetCol = 0;
            offsetRow = offsetRow + m.getNumRows();
        }

        SCAI_LOG_DEBUG( logger, "offsets for next disassembling: row " << offsetRow << ", col " << offsetCol)
    }

    auto repColDist = std::make_shared<NoDistribution>( colDist->getGlobalSize() );

    allocate( rowDist, colDist ); 
    
    fillFromAssembly( assembly );

    redistribute( rowDist, colDist );  // apply column distribution for halo computation
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void Matrix<ValueType>::vcat( const Matrix<ValueType>& m1, const Matrix<ValueType>& m2 )
{
    SCAI_ASSERT_EQ_ERROR( m1.getRowDistribution(), m2.getRowDistribution(), "vcat: matrices must have same row distribution" )

    DistributionPtr rowDist = m1.getRowDistributionPtr();

    DistributionPtr colDist( new NoDistribution( m1.getNumColumns() + m2.getNumColumns() ) );

    std::vector<const Matrix<ValueType>*> matrices;

    matrices.push_back( &m1 );
    matrices.push_back( &m2 );

    concatenate( rowDist, colDist, matrices );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void Matrix<ValueType>::hcat( const Matrix<ValueType>& m1, const Matrix<ValueType>& m2 )
{
    SCAI_ASSERT_EQ_ERROR( m1.getNumColumns(), m2.getNumColumns(), "No horizontal cut possible due to different column sizes" )

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr rowDist( new BlockDistribution( m1.getNumRows() + m2.getNumRows(), comm ) );
    DistributionPtr colDist( new NoDistribution( m1.getNumColumns() ) );

    std::vector<const Matrix<ValueType>*> matrices;

    matrices.push_back( &m1 );
    matrices.push_back( &m2 );

    concatenate( rowDist, colDist, matrices );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void Matrix<ValueType>::setRawDenseData( const IndexType n, const IndexType m, const ValueType* values )
{
    // Note: by using HArrayRef a copy of the input data is completely avoided.

    DenseStorage<ValueType> denseStorage( n, m, hmemo::HArrayRef<ValueType>( n * m, values ) );

    // virtual method makes sure that each matrix class uses his own efficient way of assignment conversion

    assign( denseStorage );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::fillFromAssembly( const MatrixAssembly<ValueType>& assembly, common::BinaryOp op )
{
    DistributionPtr colDist = getColDistributionPtr();
    DistributionPtr rowDist = getRowDistributionPtr();

    IndexType numColumns = colDist->getGlobalSize();

    if ( !colDist->isReplicated() )
    {
        // join local/halo 

        redistribute( rowDist, std::make_shared<NoDistribution>( numColumns ) );
    }

    // Note: local storage contains all owned elements

    MatrixStorage<ValueType>& localStorage = const_cast<MatrixStorage<ValueType>&>( getLocalStorage() );

    COOStorage<ValueType> cooLocal = assembly.buildLocalCOO( *rowDist, numColumns );

    hmemo::HArray<IndexType> cooIA;
    hmemo::HArray<IndexType> cooJA;
    hmemo::HArray<ValueType> cooValues;

    IndexType dummyM;
    IndexType dummyN;

    cooLocal.splitUp( dummyM, dummyN, cooIA, cooJA, cooValues );

    localStorage.fillCOO( std::move( cooIA ), std::move( cooJA ), std::move( cooValues ), op );

    if ( !colDist->isReplicated() )
    {
        // split local/halo 

        redistribute( rowDist, colDist );
    }
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Matrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */

