/**
 * @file DenseUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Implementation and instantiation of utility methods for dense storage.
 * @author Thomas Brandes
 * @date 29.05.2018
 */

#include <scai/sparsekernel/DenseUtils.hpp>

#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/Constants.hpp>

namespace scai
{

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( DenseUtils::logger, "DenseUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::convertDense2CSR(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    getSparseRowSizes( csrIA, numRows, numColumns, denseValues, prefLoc );

    SCAI_REGION( "Sparse.Dense.buildCSR" )

    IndexType numValues = HArrayUtils::scan1( csrIA, prefLoc );

    static LAMAKernel<DenseKernelTrait::getCSRValues<ValueType> > getCSRValues;

    ContextPtr loc = prefLoc;
    getCSRValues.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    WriteOnlyAccess<IndexType> wJA( csrJA, loc, numValues );
    WriteOnlyAccess<ValueType> wValues( csrValues, loc, numValues );

    ReadAccess<ValueType> rDense( denseValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    RealType<ValueType> eps = 0;

    getCSRValues[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns,
                       rDense.get(), eps );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::getSparseRowSizes(
    HArray<IndexType>& csrSizes,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( denseValues.size(), numRows * numColumns, 
                          "denseValues: size not " << numRows << " x " << numColumns )

    static LAMAKernel<DenseKernelTrait::getCSRSizes<ValueType> > getCSRSizes;

    ContextPtr loc = prefLoc;
    getCSRSizes.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rDenseValues( denseValues, loc );

    // allocate sizes with one more element, so it can be used for offsets
    WriteOnlyAccess<IndexType> wSizes( csrSizes, loc, numRows + 1 );

    RealType<ValueType> eps = 0;

    getCSRSizes[loc]( wSizes.get(), numRows, numColumns, rDenseValues.get(), eps );
    wSizes.resize( numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::convertCSR2Dense(
    HArray<ValueType>& denseValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.Dense.setCSR" )

    static LAMAKernel<DenseKernelTrait::setCSRValues<ValueType> > setCSRValues;

    ContextPtr loc = prefLoc;
    setCSRValues.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );

    WriteOnlyAccess<ValueType> wDense( denseValues, loc, numRows * numColumns );

    setCSRValues[loc]( wDense.get(), numRows, numColumns, rIA.get(), rJA.get(), rValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const IndexType n,
    const HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( n, oldSolution.size(), "size mismatch" )
    SCAI_ASSERT_EQ_DEBUG( n, rhs.size(), "size mismatch" )

    static LAMAKernel<DenseKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rValues( denseValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, n );
    jacobi[loc]( wSolution.get(), n, rValues.get(),
                 rOldSolution.get(), rRhs.get(), omega );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::jacobiHalo(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& diagonal,
    const HArray<ValueType>& oldSolution,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    static LAMAKernel<DenseKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = prefLoc;
    jacobiHalo.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rValues( denseValues, loc ); 
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numRows );

    jacobiHalo[loc]( wSolution.get(), rDiagonal.get(),
                     numRows, numColumns, rValues.get(),
                     rOldSolution.get(), omega );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* DenseUtils::gemv0(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& denseValues,
    const common::MatrixOp op,
    const bool,
    ContextPtr prefLoc )
{
    // Initialize result array in any case with 0

    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), prefLoc );

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        return NULL;   // already done
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;

    ContextPtr loc = prefLoc;
    gemv.getSupportedContext( loc );

    const IndexType lda = numColumns; // stride for denseValues between rows
    const IndexType ldx = 1;
    const IndexType ldresult = 1;

    ReadAccess<ValueType> rDenseValues( denseValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );

    SCAI_CONTEXT_ACCESS( loc )

    gemv[loc]( op, numRows, numColumns, alpha, rDenseValues.get(), lda, 
               rX.get(), ldx, ValueType( 0 ), wResult.get(), ldresult );

    return NULL;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* DenseUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& denseValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // if beta is 0, call the simpler routine, avoids accesses to y

    if ( beta == common::Constants::ZERO )
    {
       return gemv0( result, alpha, x, numRows, numColumns, denseValues, op, async, prefLoc );
    }

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        // result = beta * y, beta != 0

        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, prefLoc );
        return NULL;
    }

    // split it up: result = y; result = alpha * denseMatrix * x + beta * result

    utilskernel::HArrayUtils::assign( result, y, prefLoc );

    // use BLAS-2 routine for result = alpha * denseMatrix * x + beta * result

    static LAMAKernel<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;

    ContextPtr loc = prefLoc;

    gemv.getSupportedContext( loc );

    const IndexType lda = numColumns; // stride for denseValues between rows
    const IndexType ldx = 1;
    const IndexType ldresult = 1;

    ReadAccess<ValueType> rDenseValues( denseValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );

    SCAI_CONTEXT_ACCESS( loc )

    gemv[loc]( op, numRows, numColumns, alpha, rDenseValues.get(), lda, 
               rX.get(), ldx, beta, wResult.get(), ldresult );

    return NULL;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::gemm(
    HArray<ValueType>& c,
    const ValueType alpha,
    const HArray<ValueType>& a,
    const common::MatrixOp opA,
    const HArray<ValueType>& b,
    const common::MatrixOp opB,
    const ValueType beta,
    const IndexType m, 
    const IndexType n, 
    const IndexType k,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( c.size(), m * n, "c has not size " << m << " x " << n )
    SCAI_ASSERT_EQ_DEBUG( a.size(), m * k, "a has not size " << m << " x " << k )
    SCAI_ASSERT_EQ_DEBUG( b.size(), k * n, "b has not size " << k << " x " << n )

    if ( m == 0 || n == 0 || k == 0 )
    {
        return;  // done for trivial cases
    }

    // A is m x k if normal or A is k x m if transposed

    const IndexType lda = common::isTranspose( opA ) ? m : k; 
    const IndexType ldb = common::isTranspose( opB ) ? k : n;  // B is k x n
    const IndexType ldc = n;  // C is m x n

    static LAMAKernel<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;

    ContextPtr loc = prefLoc;

    gemm.getSupportedContext( loc );

    ReadAccess<ValueType> aAccess( a, loc );
    ReadAccess<ValueType> bAccess( b, loc );
    WriteAccess<ValueType> resAccess( c, loc );

    SCAI_CONTEXT_ACCESS( loc )

    gemm[loc]( opA, opB,
               m, n, k, alpha, aAccess.get(), lda, bAccess.get(), ldb, beta,
               resAccess.get(), ldc );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::invert(
    HArray<ValueType>& a,
    const IndexType n,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( a.size(), n * n, "matrix a has not size " << n << " x " << n )

    static LAMAKernel<blaskernel::BLASKernelTrait::getinv<ValueType> > getinv;

    ContextPtr loc = prefLoc;
    getinv.getSupportedContext( loc );
    WriteAccess<ValueType> wA( a, loc );
    SCAI_CONTEXT_ACCESS( loc );
    const IndexType lda = n;
    getinv[loc]( n, wA.get(), lda );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::setScalar(
    HArray<ValueType>& values,
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType scalar,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    static LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;
    ContextPtr loc = prefLoc;
    setValue.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wData( values, loc );
    setValue[loc]( wData.get(), numRows, numColumns, scalar, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::setRows(
    HArray<ValueType>& denseValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& rowValues,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( rowValues.size(), numRows, "illegal size" )

    static LAMAKernel<DenseKernelTrait::setRows<ValueType> > setRows;

    ContextPtr loc = prefLoc;
    setRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rRows( rowValues, loc );
    WriteAccess<ValueType> wValues( denseValues, loc );
    setRows[loc]( wValues.get(), numRows, numColumns, rRows.get(), op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::setColumns(
    HArray<ValueType>& denseValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<ValueType>& columnValues,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( columnValues.size(), numColumns, "illegal size" )

    static LAMAKernel<DenseKernelTrait::setColumns<ValueType> > setColumns;

    ContextPtr loc = prefLoc;
    setColumns.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rColumns( columnValues, loc );
    WriteAccess<ValueType> wValues( denseValues, loc );
    setColumns[loc]( wValues.get(), numRows, numColumns, rColumns.get(), op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DenseUtils::getNumValues(
    const HArray<ValueType>& denseValues,
    const IndexType numRows,
    const IndexType numColumns,
    ContextPtr prefLoc )
{
    static LAMAKernel<DenseKernelTrait::nonZeroValues<ValueType> > nonZeroValues;

    ContextPtr loc = prefLoc;
    nonZeroValues.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rValues( denseValues, loc );
    ValueType zero = 0;
    IndexType count = nonZeroValues[loc]( rValues.get(), numRows, numColumns, zero );
    return count;
}

/* -------------------------------------------------------------------------- */

#define DENSE_UTILS_SPECIFIER( ValueType )           \
                                                     \
    template void DenseUtils::getSparseRowSizes(     \
        HArray<IndexType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::convertDense2CSR(      \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::convertCSR2Dense(      \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::jacobi(                \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::jacobiHalo(            \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template tasking::SyncToken* DenseUtils::gemv(   \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        const common::MatrixOp,                      \
        const bool,                                  \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::gemm(                  \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const common::MatrixOp,                      \
        const HArray<ValueType>&,                    \
        const common::MatrixOp,                      \
        const ValueType,                             \
        const IndexType,                             \
        const IndexType,                             \
        const IndexType,                             \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::invert(                \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        ContextPtr );                                \
                                                     \
    template IndexType DenseUtils::getNumValues(     \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::setScalar(             \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const ValueType,                             \
        const common::BinaryOp op,                   \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::setRows(               \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        const common::BinaryOp op,                   \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::setColumns(            \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        const common::BinaryOp op,                   \
        ContextPtr );                                \

SCAI_COMMON_LOOP( DENSE_UTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef DENSE_UTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
