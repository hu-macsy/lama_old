/**
 * @file DenseStorage.cpp
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
 * @brief Instantiation for template class DenseStorage.
 * @author Thomas Brandes, Michael Drost
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/DenseUtils.hpp>
#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/freeFunction.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>
#include <cmath>


using std::abs;
// so we can use abs for float and double and own abs for Complex

using std::shared_ptr;

namespace scai
{

using common::TypeTraits;

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using utilskernel::UtilKernelTrait;

using sparsekernel::DenseKernelTrait;
using sparsekernel::CSRUtils;
using sparsekernel::DenseUtils;

using common::BinaryOp;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseStorage<ValueType>::logger,
                              "MatrixStorage.DenseStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
HArray<ValueType>& DenseStorage<ValueType>::getData()
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& DenseStorage<ValueType>::getValues() const
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DenseStorage<ValueType>::getNumValues() const
{
    static LAMAKernel<DenseKernelTrait::nonZeroValues<ValueType> > nonZeroValues;
    ContextPtr loc = getContextPtr();
    nonZeroValues.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> values( mData, loc );
    ValueType zero = 0;
    IndexType count = nonZeroValues[loc]( values.get(), getNumRows(), getNumColumns(), zero );
    SCAI_LOG_INFO( logger, *this << ": #non-zero values = " << count )
    return count;
}

/* --------------------------------------------------------------------------- */

#ifdef SCAI_ASSERT_LEVEL_OFF
template<typename ValueType>
void DenseStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void DenseStorage<ValueType>::check( const char* /* msg */ ) const
{
    SCAI_ASSERT_EQUAL_ERROR( getNumRows() * getNumColumns(), mData.size() )
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    // ToDo: avoid temporary array row by new version buildSparseArray with offs and inc argument

    HArray<ValueType> row;
    getRow( row, i );
    ValueType zero = 0;
    HArrayUtils::buildSparseArray( values, jA, row, zero, mContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const
{
    HArray<ValueType> col;
    getColumn( col, j );
    ValueType zero = 0;
    HArrayUtils::buildSparseArray( values, iA, col, zero, mContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType rowIndex ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rowIndex, getNumRows(), "row index out of range" )

    row.clear();                    // make all data invalid
    row.resize( getNumColumns() );      // resize it

    // inc = denseindex( i, j + 1, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( rowIndex, 0, numRows, numColumns )

    const IndexType inc   = 1;
    const IndexType first = rowIndex * getNumColumns();

    HArrayUtils::setArraySection( row, 0, 1,             // row (:)
                                  mData, first, inc,
                                  getNumColumns(), BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType rowIndex,
        const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rowIndex, getNumRows(), "row index out of range" )

    // inc = denseindex( i, j + 1, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( i, 0, numRows, numColumns )

    const IndexType inc   = 1;
    const IndexType first = rowIndex * getNumColumns();

    HArrayUtils::setArraySection( mData, first, inc,
                                  row, 0, 1,
                                  getNumColumns(),
                                  op,
                                  getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    column.clear();                 // make all data invalid
    column.resize( getNumRows() );      // resize it

    // inc = denseindex( i + 1, j, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( 0, j, numRows, numColumns )

    const IndexType inc   = getNumColumns();
    const IndexType first = j;

    HArrayUtils::setArraySection( column, 0, 1,
                                  mData, first, inc,
                                  getNumRows(),
                                  BinaryOp::COPY,
                                  getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    // inc = denseindex( i + 1, j, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( 0, j, numRows, numColumns )

    const IndexType inc   = getNumColumns();
    const IndexType first = j;

    HArrayUtils::setArraySection( mData, first, inc ,
                                  column, 0, 1,
                                  getNumRows(),
                                  op,
                                  getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    IndexType numDiagonalValues = common::Math::min( getNumColumns(), getNumRows() );

    diagonal.clear();                       // make all data invalid
    diagonal.resize( numDiagonalValues );   // resize it

    // inc = denseindex( i + 1, i + 1, numRows, numColumns ) - denseindex( i, i, numRows, numColumns )
    // first = denseindex( 0, 0, numRows, numColumns )

    const IndexType inc   = getNumRows() + 1;
    const IndexType first = 0;

    HArrayUtils::setArraySection( diagonal, 0, 1,
                                  mData, first, inc,
                                  numDiagonalValues,
                                  BinaryOp::COPY,
                                  mData.getValidContext() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setDiagonal( const ValueType value )
{
    const IndexType numDiagonalValues = common::Math::min( getNumColumns(), getNumRows() );

    // inc = denseindex( i + 1, i + 1, numRows, numColumns ) - denseindex( i, i, numRows, numColumns )
    // first = denseindex( 0, 0, numRows, numColumns )
    // dense values are stored row-wise

    const IndexType inc   = getNumRows() + 1;
    const IndexType first = 0;

    HArrayUtils::fillArraySection( mData, first, inc,
                                   value, numDiagonalValues,
                                   BinaryOp::COPY,
                                   mData.getValidContext() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    const IndexType numDiagonalValues = common::Math::min( getNumColumns(), getNumRows() );

    SCAI_ASSERT_GE_DEBUG( diagonal.size(), numDiagonalValues, "diagonal array has insufficient size" )

    // inc = denseindex( i + 1, i + 1, numRows, numColumns ) - denseindex( i, i, numRows, numColumns )
    // first = denseindex( 0, 0, numRows, numColumns )

    const IndexType inc   = getNumRows() + 1;
    const IndexType first = 0;

    HArrayUtils::setArraySection( mData, first, inc,
                                  diagonal, 0, 1,
                                  numDiagonalValues,
                                  BinaryOp::COPY,
                                  mData.getValidContext() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::scale( const ValueType value )
{
    // not used here: HArrayUtils::scale( mData, value, this->getContextPtr() )
    // reasoning:     workload distribution would not fit to distribution of rows
    static LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;
    ContextPtr loc = getContextPtr();
    setValue.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wData( mData, loc );
    setValue[loc]( wData.get(), getNumRows(), getNumColumns(), value, BinaryOp::MULT );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mData, mData, common::UnaryOp::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::scaleRows( const HArray<ValueType>& values )
{
    static LAMAKernel<DenseKernelTrait::scaleRows<ValueType> > denseScaleRows;
    ContextPtr loc = this->getContextPtr();
    denseScaleRows.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rDiagonal( values, loc );
    WriteAccess<ValueType> wData( mData, loc );
    denseScaleRows[loc]( wData.get(), getNumRows(), getNumColumns(), rDiagonal.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::transposeImpl()
{
    SCAI_REGION( "Storage.Dense.transpose" )

    HArrayUtils::transpose( mData, getNumColumns(), getNumRows(), mData, false, getContextPtr() );
 
    _MatrixStorage::setDimension( getNumColumns(), getNumRows() );
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DenseStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( ValueType ) * mData.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType DenseStorage<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format DenseStorage<ValueType>::getFormat() const
{
    return Format::DENSE;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setIdentity( const IndexType size )
{
    allocate( size, size );            // sets also all values to zero
    setDiagonal( ValueType( 1 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    allocate( size, size );     // sets also all values to zero
    setDiagonalV( diagonal );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setZero()
{
    SCAI_ASSERT_EQ_DEBUG( mData.size(), getNumRows() * getNumColumns(), 
                          "illegal data size for DenseStorage " << getNumRows() << " x " << getNumColumns() )
    ValueType zero = 0;

    HArrayUtils::setScalar( mData, zero, common::BinaryOp::COPY, this->getContextPtr() );

    SCAI_LOG_DEBUG( logger, *this << " has been set to zero" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
{
    DenseUtils::getSparseRowSizes( ia, getNumRows(), getNumColumns(), mData, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::buildCSRData( 
    hmemo::HArray<IndexType>& csrIA, 
    hmemo::HArray<IndexType>& csrJA, 
    hmemo::_HArray& csrValues ) const
{
    if ( csrValues.getValueType() == getValueType() )
    {
        HArray<ValueType>& typedCSRValues = static_cast<HArray<ValueType>&>( csrValues );

        DenseUtils::convertDense2CSR( csrIA, csrJA, typedCSRValues,
                                      getNumRows(), getNumColumns(), mData, getContextPtr() );
    }
    else
    {
        HArray<ValueType> tmpValues;

        DenseUtils::convertDense2CSR( csrIA, csrJA, tmpValues, 
                                      getNumRows(), getNumColumns(), mData, getContextPtr() );

        HArrayUtils::_assign( csrValues, tmpValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setCSRData(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    if ( values.getValueType() == getValueType() )
    {
        setCSRDataImpl( numRows, numColumns, ia, ja, 
		        static_cast<const HArray<ValueType>&>( values ) );
    }
    else
    {
        setCSRDataImpl( numRows, numColumns, ia, ja, 
                        utilskernel::convertHArray<ValueType>( values, getContextPtr() ) );
    }
}

template<typename ValueType>
void DenseStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )
{   
    IndexType numValues = ja.size();

    
    SCAI_LOG_INFO( logger,
                   "setCRSData for dense storage " << numRows << " x " << numColumns << ", nnz = " << numValues )
    
    _MatrixStorage::setDimension( numRows, numColumns );
    
    if ( ia.size() == numRows )
    {   
        HArray<IndexType> tmpOffsets;
        IndexType total = CSRUtils::sizes2offsets( tmpOffsets, ia, getContextPtr() );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        setCSRDataImpl( numRows, numColumns, tmpOffsets, ja, values );
        return;
    }
    
    SCAI_ASSERT_EQ_ERROR( ia.size(), numRows + 1, "size mismatch of csr IA array" )
    SCAI_ASSERT_ERROR( CSRUtils::validOffsets( ia, numValues, getContextPtr() ), "illegal CSR offset array" );
    
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( ja, numColumns, getContextPtr() ),
                       "CSR ja array contains illegal column indexes, #columns = " << numColumns );
    
    DenseUtils::convertCSR2Dense( mData, numRows, numColumns, ia, ja, values, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invert( " << other << ") to a dense storage" )

    if ( other.getFormat() == Format::DENSE )
    {
        const DenseStorage<ValueType>* otherDense = dynamic_cast<const DenseStorage<ValueType>*>( &other );
        SCAI_ASSERT_ERROR( otherDense, "Internal error: dynamic cast Dense" )
        invertDense( *otherDense );
    }
    else
    {
        SCAI_UNSUPPORTED( "invert (" << other << ") requires conversion to dense storage" )
        assignImpl( other );
        invertDense( *this ); // is always done in place
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::invertDense( const DenseStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invertDense: " << other )

    // invert is always done in place, so assign other to this

    if ( &other != this )
    {
        SCAI_LOG_INFO( logger, "invertDense: copy input matrix to this matrix" )
        assignDense( other );
    }

    int nRows = this->getNumRows();
    int nCols = this->getNumColumns();
    SCAI_ASSERT_EQUAL_ERROR( nRows, nCols )
    static LAMAKernel<blaskernel::BLASKernelTrait::getinv<ValueType> > getinv;
    ContextPtr loc = this->getContextPtr();
    getinv.getSupportedContext( loc );
    WriteAccess<ValueType> denseValues( mData, loc );
    SCAI_CONTEXT_ACCESS( loc );
    getinv[loc]( nRows, denseValues.get(), nCols );
    SCAI_LOG_INFO( logger, "invertDense: this = " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::binaryOp(
    const MatrixStorage<ValueType>& a,
    const common::BinaryOp op,
    const MatrixStorage<ValueType>& b )
{
    if ( a.getFormat() != Format::DENSE )
    {
        if ( &b == this )
        {
            auto aDense = convert<DenseStorage<ValueType>>( a );
            binaryOp( aDense, op, b );
        }
        else
        {
            // reuse this storage for conversion of a
            assign( a );
            binaryOp( *this, op, b );
        }
    }
    else if ( b.getFormat() != Format::DENSE )
    {
        if ( &a == this )
        {
            auto bDense = convert<DenseStorage<ValueType>>( b );
            binaryOp( a, op, bDense );
        }
        else
        {
            // reuse this storage for conversion of a
            assign( b );
            binaryOp( a, op, *this );
        }
    }
    else
    {
        binaryOpDense( static_cast<const DenseStorage<ValueType>&>( a ), op,
                       static_cast<const DenseStorage<ValueType>&>( b ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::binaryOpDense(
    const DenseStorage<ValueType>& a,
    const common::BinaryOp op, 
    const DenseStorage<ValueType>& b )
{
    SCAI_ASSERT_EQ_ERROR( a.getNumRows(), b.getNumRows(), "serious size mismatch" )
    SCAI_ASSERT_EQ_ERROR( a.getNumColumns(), b.getNumColumns(), "serious size mismatch" )

    _MatrixStorage::setDimension( a.getNumRows(), a.getNumColumns() );

    utilskernel::HArrayUtils::binaryOp( mData, a.getValues(), b.getValues(), op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    const IndexType nSource = common::isTranspose( op ) ? getNumRows() : getNumColumns();
    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    if ( alpha == common::Constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously

        if ( beta == common::Constants::ZERO )
        {
            HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), this->getContextPtr() );
        }
        else
        {
            HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
        }

        return;
    }

    SCAI_ASSERT_EQ_ERROR( x.size(), nSource, "vector x in A * x has illegal size" )

    if ( beta != common::Constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "size mismatch y, beta = " << beta )
    }

    if ( nTarget == 0 )
    {
        result.clear();  // result will get size zero
        return;          // nothing to do
    }

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector, try on " << *mContext )

    // using BLAS2 interface requires result and y to be aliased

    if ( beta == common::Constants::ZERO )
    {
        utilskernel::HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), getContextPtr() );
    }
    else if ( &result != &y )
    {
        utilskernel::HArrayUtils::assign( result, y, this->getContextPtr() );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if ( nSource == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if ( beta == common::Constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if ( beta == common::Constants::ONE )
        {
            // no scaling required
        }
        else
        {
            HArrayUtils::compute( result, beta, BinaryOp::MULT, result, this->getContextPtr() );
        }
    }
    else
    {
        // mNumColums > 0, mnumRows > 0, so we avoid problems for gemv with m==0 or n==0
        static LAMAKernel<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
        ContextPtr loc = this->getContextPtr();
        gemv.getSupportedContext( loc );
        ReadAccess<ValueType> denseValues( mData, loc );
        ReadAccess<ValueType> rX( x, loc );
        int lda = getNumColumns(); // stride for denseValues between rows
        WriteAccess<ValueType> wResult( result, loc );
        SCAI_CONTEXT_ACCESS( loc )
        // gemv:  result = alpha * this * x + beta * result
        gemv[loc]( op, getNumRows(), getNumColumns(), alpha, denseValues.get(), lda, rX.get(),
                   1, beta, wResult.get(), 1 );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    ReadAccess<ValueType> values( mData );
    stream << "DenseStorage " << getNumRows() << " x " << getNumColumns() << ", addr  = " << values.get() << endl;

    for ( IndexType i = 0; i < getNumRows(); i++ )
    {
        stream << "Row " << i << " :";

        for ( IndexType j = 0; j < getNumColumns(); j++ )
        {
            stream << " " << values[i * getNumColumns() + j];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    SCAI_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B + " << beta << " * C" << ", with A = " << a << ", B = " << b << ", C = " << c )

    if ( a.getFormat() != Format::DENSE )
    {
        matrixTimesMatrix( alpha, convert<DenseStorage<ValueType>>( a ), b, beta, c );
    }
    else if ( b.getFormat() != Format::DENSE )
    {
        matrixTimesMatrix( alpha, a, convert<DenseStorage<ValueType>>( b ), beta, c );
    }
    else if ( c.getFormat() != Format::DENSE )
    {
        if ( beta == 0 )
        {
            matrixTimesMatrix( alpha, a, b, beta, zero<DenseStorage<ValueType>>( a.getNumRows(), b.getNumColumns() )  );
        }
        else
        {
            matrixTimesMatrix( alpha, a, b, beta, convert<DenseStorage<ValueType>>( c )  );
        }
    }
    else
    {
        // all matrices are now dense
        matrixTimesMatrixDense( alpha, static_cast<const DenseStorage<ValueType>&>( a ), 
                                       static_cast<const DenseStorage<ValueType>&>( b ),
                                beta,  static_cast<const DenseStorage<ValueType>&>( c ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesMatrixDense(
    const ValueType alpha,
    const DenseStorage<ValueType>& a,
    const DenseStorage<ValueType>& b,
    const ValueType beta,
    const DenseStorage<ValueType>& c )
{
    // shape(a) = m x k,  shape(b) = k x n
    const DenseStorage<ValueType>* ptrA = &a;
    const DenseStorage<ValueType>* ptrB = &b;
    std::shared_ptr<DenseStorage<ValueType> > tmpA;
    std::shared_ptr<DenseStorage<ValueType> > tmpB;
    SCAI_LOG_INFO( logger,
                   "matrixTimesMatrixDense: " << alpha << " * a * b + " << beta << " * c, with a = " << a << ", b = " << b << ", c = " << c )

    if ( &a == this )
    {
        SCAI_LOG_INFO( logger, "temporary for A in A * B ( dense storages) needed" )
        tmpA = std::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( a ) );
        ptrA = tmpA.get();
    }

    if ( &b == this )
    {
        // avoid two temporaries
        if ( &a == &b )
        {
            SCAI_LOG_INFO( logger, "same temporary for A and B in A * B ( dense storages) needed" )
            ptrB = tmpA.get();
        }
        else
        {
            SCAI_LOG_INFO( logger, "temporary for B in A * B ( dense storages) needed" )
            tmpB = std::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( b ) );
            ptrB = tmpB.get();
        }
    }

    IndexType m = a.getNumRows();
    IndexType k = b.getNumRows();
    IndexType n = b.getNumColumns();

    SCAI_ASSERT_EQUAL_ERROR( k, a.getNumColumns() )
 
    _MatrixStorage::setDimension( m, n );

    if ( beta == common::Constants::ZERO )
    {
        // do not care at all about C as it might be any dummy, or aliased to result
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );
        SCAI_LOG_INFO( logger, "init this result with 0, size = " << m * n )
        WriteOnlyAccess<ValueType> resAccess( mData, loc, m * n );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( resAccess.get(), m * n, ValueType( 0 ), BinaryOp::COPY );
    }
    else if ( this != &c )
    {
        // force result = C
        SCAI_ASSERT_EQUAL_ERROR( m, c.getNumRows() )
        SCAI_ASSERT_EQUAL_ERROR( n, c.getNumColumns() )

        _MatrixStorage::setDimension( m, n );

        static LAMAKernel<blaskernel::BLASKernelTrait::copy<ValueType> > copy;
        ContextPtr loc = this->getContextPtr();
        copy.getSupportedContext( loc );
        ReadAccess<ValueType> cAccess( c.getValues(), loc );
        WriteOnlyAccess<ValueType> resAccess( mData, loc, m * n );
        SCAI_LOG_TRACE( logger, "Copying: res = c " )
        SCAI_CONTEXT_ACCESS( loc )
        copy[loc]( n * m, cAccess.get(), 1, resAccess.get(), 1 );
    }
    else
    {
        SCAI_LOG_INFO( logger, "results is aliased with C as required for gemm, beta = " << beta )
    }

    // now C used for BLAS3 call of GEMM is this matrix
    int lda = a.getNumColumns();
    int ldb = b.getNumColumns();
    int ldc = getNumColumns();
    // Now  C = alpha * A * B + beta * C, can use GEMM of BLAS3
    SCAI_LOG_TRACE( logger,
                    "GEMM( Normal, Normal, m:" << m << ", n:" << n << ", k:" << k << ", alpha:" << alpha << ", A , lda:" << lda << " ,B ,ldb:" << ldb << " ,beta:" << beta << " C ,ldc:" << ldc << " )" )

    if ( lda != 0 && n != 0 && m != 0 )
    {
        static LAMAKernel<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;
        ContextPtr loc = this->getContextPtr();
        gemm.getSupportedContext( loc );
        ReadAccess<ValueType> aAccess( ptrA->getValues(), loc );
        ReadAccess<ValueType> bAccess( ptrB->getValues(), loc );
        WriteAccess<ValueType> resAccess( mData, loc );
        SCAI_CONTEXT_ACCESS( loc )
        gemm[loc]( common::MatrixOp::NORMAL, common::MatrixOp::NORMAL, 
                   m, n, k, alpha, aAccess.get(), lda, bAccess.get(), ldb, beta,
                   resAccess.get(), ldc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::~DenseStorage()
{
    SCAI_LOG_DEBUG( logger, "~DenseStorage for matrix " << getNumRows() << " x " << getNumColumns() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseStorage<ValueType>::l1Norm() const
{
    IndexType n = getNumRows() * getNumColumns();

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    ContextPtr loc = this->getContextPtr();
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> data( mData, loc );
    SCAI_CONTEXT_ACCESS( loc );
    IndexType inc = 1;
    return asum[loc]( n, data.get(), inc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseStorage<ValueType>::l2Norm() const
{
    IndexType n = getNumRows() * getNumColumns();

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = this->getContextPtr();
    dot.getSupportedContext( loc );
    ReadAccess<ValueType> data( mData, loc );
    SCAI_CONTEXT_ACCESS( loc );
    return common::Math::sqrt( dot[loc]( n, data.get(), 1, data.get(), 1 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseStorage<ValueType>::maxNorm() const
{
    IndexType n = getNumRows() * getNumColumns();

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    return HArrayUtils::maxNorm( mData, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_REGION( "Storage.Dense.maxDiffNorm" )

    SCAI_ASSERT_EQ_ERROR( getNumRows(), other.getNumRows(), "row size mismatch for maxDiffNorm" )
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), other.getNumColumns(), "col size mismatch for maxDiffNorm" )

    SCAI_LOG_INFO( logger, *this << ": maxDiffNorm( " << other << " )" )

    if ( other.getFormat() == Format::DENSE )
    {
        return maxDiffNormImpl( static_cast<const DenseStorage<ValueType>&>( other ) );
    }
    else
    {
        SCAI_UNSUPPORTED( other << ": converted to " << typeName() << " for maxDiffNorm" )

        return maxDiffNormImpl( convert<DenseStorage<ValueType>>( other ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseStorage<ValueType>::maxDiffNormImpl( const DenseStorage<ValueType>& other ) const
{
    // no more checks needed here

    IndexType n = getNumRows() * getNumColumns();

    if ( n == 0 )
    {
        return static_cast<RealType<ValueType> >( 0 );
    }

    return HArrayUtils::maxDiffNorm( mData, other.mData, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::assignDense( const DenseStorage<OtherValueType>& other )
{
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {
        SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )
        return;
    }

    // actualize member variables of base class

    _MatrixStorage::_assign( other ); // copy sizes, flags

    // copy the dense data, maybe with type conversion

    HArrayUtils::setArray( mData, other.getValues(), common::BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<DenseStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    if ( other.getFormat() == Format::DENSE )
    {
        // both storage have COO format, use special method for it

        assignDense( static_cast<const DenseStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRData( otherCSR.getNumRows(), otherCSR.getNumColumns(),
                    otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues() );

        SCAI_LOG_INFO( logger, "DenseStorage: assign CSR : " << other << ", this = " << *this );
    }
    else
    {
        HArray<IndexType>  csrIA;
        HArray<IndexType>  csrJA;
        HArray<ValueType>  csrValues;     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        // just a thought for optimization: use mIA, mJA, mValues instead of csrIA, csrJA, csrValues
        // but does not help much at all as resort of entries requires already temporaries.

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );

        SCAI_LOG_INFO( logger, "DenseStorage: assign " << other << ", built CSR tmp data with ia = " 
                                << csrIA << ", ja = " << csrJA << ", values = " << csrValues << ", this = " << *this )

    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate dense storage of size " << numRows << " x " << numColumns )
  
    _MatrixStorage::setDimension( numRows, numColumns );

    mData.resize( getNumRows() * getNumColumns() );

    setZero();

    SCAI_LOG_DEBUG( logger, *this << " allocated, #values = " << getNumRows() * getNumColumns() << ", not initialized" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << getTypeName() << "( rows = " << getNumRows() << ", cols = " << getNumColumns() << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    const IndexType pos = i * getNumColumns() + j;

    return utilskernel::HArrayUtils::getVal<ValueType>( mData, pos );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setValue( const IndexType i,
                                        const IndexType j,
                                        const ValueType val,
                                        const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    const IndexType pos = i * getNumColumns() + j;

    utilskernel::HArrayUtils::setVal( mData, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mData.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::wait() const
{
    mData.wait();
}

/* --------------------------------------------------------------------------- */
/*  Constructors for DenseStorage                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( ContextPtr ctx ) :

    MatrixStorage<ValueType>( 0, 0, ctx ),
    mData( ctx )
{
}

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const IndexType numRows, const IndexType numColumns, ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx )

{
    mData.resize( numRows * numColumns );
    this->setZero();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage(
    const IndexType numRows,
    const IndexType numColumns,
    HArray<ValueType> denseData, 
    ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mData( std::move( denseData ) )

{
    SCAI_ASSERT_EQ_ERROR( mData.size(), numRows * numColumns, 
                          "size of array does not fit the shape " << numRows << " x " << numColumns )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const DenseStorage<ValueType>& other ) :

    MatrixStorage<ValueType>( other )

{
    assignDense( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( DenseStorage<ValueType>&& other ) :

   MatrixStorage<ValueType>( std::move( other ) ),
   mData( std::move( other.mData ) )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>& DenseStorage<ValueType>::operator=( const DenseStorage<ValueType>& other )
{
    assignDense( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>& DenseStorage<ValueType>::operator=( DenseStorage<ValueType>&& other )
{
    // move of all member variables
    
    mData = std::move( other.mData );

    // call of move assignment for base class 

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    // Note: other ends up in a zero storage

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::swap( DenseStorage<ValueType>& other )
{
    MatrixStorage<ValueType>::swap( other ); // swap member variable of base class

    mData.swap( other.mData );               // swap data
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::splitUp( IndexType& numRows, IndexType& numColumns, HArray<ValueType>& values )
{
    values = std::move( mData );
    _MatrixStorage::splitUp( numRows, numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>*
DenseStorage<ValueType>::copy() const
{
    return new DenseStorage<ValueType>( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
DenseStorage<ValueType>* DenseStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<DenseStorage<ValueType> > storage( new DenseStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

template<typename ValueType>
std::string DenseStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "DenseStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* DenseStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* DenseStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType DenseStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DENSE, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType DenseStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* DenseStorage<ValueType>::create()
{
    return new DenseStorage<ValueType>();
}

/* ========================================================================= */

template<typename ValueType>
void DenseStorage<ValueType>::fillCOO(
    hmemo::HArray<IndexType> ia,
    hmemo::HArray<IndexType> ja,
    hmemo::HArray<ValueType> values,
    const common::BinaryOp op )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "COO data: ia and ja must have same size" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "COO data: ia and values must have same size" )

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( ia, getNumRows() ), "illegal index in ia of COO data" )
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( ja, getNumColumns() ), "illegal index in ja of COO data" )

    hmemo::HArray<IndexType> pos;  // becomes: ia * numColumns + ja, positions for dense data

    bool unique = false;    // there might be multiple entries in the COO data

    HArrayUtils::arrayPlusArray<IndexType>( pos, getNumColumns(), ia, 1, ja, getContextPtr() );
    HArrayUtils::scatter<ValueType, ValueType>( mData, pos, unique, values, op, getContextPtr() );
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
