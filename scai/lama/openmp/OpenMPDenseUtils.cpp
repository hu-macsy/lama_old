/**
 * @file OpenMPDenseUtils.cpp
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
 * @brief Implementation of Dense utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/openmp/OpenMPDenseUtils.hpp>

// local library
#include <scai/lama/DenseKernelTrait.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

namespace lama
{

using common::TypeTraits;

SCAI_LOG_DEF_LOGGER( OpenMPDenseUtils::logger, "OpenMP.DenseUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
IndexType OpenMPDenseUtils::nonZeroValues(
    const DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType eps )
{
    IndexType count = 0;

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE ) reduction( + : count )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            if ( common::Math::abs( denseValues[denseindex( i, j, numRows, numColumns )] ) > eps )
            {
                count++;
            }
        }
    }
 
    return count;
}

template<typename DenseValueType>
void OpenMPDenseUtils::getCSRSizes(
    IndexType csrSizes[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType denseValues[],
    const DenseValueType eps )
{
    if( numRows > 0 )
    {
        SCAI_ASSERT_DEBUG( csrSizes != NULL, "csrSizes is NULL" )

        if( numColumns > 0 )
        {
            SCAI_ASSERT_DEBUG( denseValues != NULL, "denseValues is NULL" )
        }
    }

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; ++i )
    {
        IndexType nonZeros = 0; // count for each row in parallel

        for( IndexType j = 0; j < numColumns; ++j )
        {
            const DenseValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            if( common::Math::abs( value ) > eps )
            {
                ++nonZeros;
            }
            else if( i == j && diagonalFlag )
            {
                ++nonZeros; // count also zero elements for diagonals
            }
        }

        csrSizes[i] = nonZeros;
    }
}

/** Helper routine for conversion Dense to CSR */

template<typename DenseValueType,typename CSRValueType>
void OpenMPDenseUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType denseValues[],
    const DenseValueType eps )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<DenseValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", size is " << numRows << " x " << numColumns )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = csrIA[i];

        if( i < numColumns && diagonalFlag )
        {
            // start with diagonal element in any case

            const DenseValueType& value = denseValues[denseindex( i, i, numRows, numColumns )];

            csrValues[offset] = static_cast<CSRValueType>( value );
            csrJA[offset] = i;
            offset++;
        }

        for( IndexType j = 0; j < numColumns; ++j )
        {
            if( i == j && diagonalFlag )
            {
                continue; // diagonal element is already first element in ja, values
            }

            const DenseValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            if( common::Math::abs( value ) > eps )
            {
                csrValues[offset] = static_cast<CSRValueType>( value );
                csrJA[offset] = j;
                offset++;
            }
        }

        // verification that offset array was really a good one
        // check is not needed if non-zero values have been counted by getCSRSizes

        SCAI_ASSERT_EQUAL_DEBUG( offset, csrIA[i + 1] )
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType,typename CSRValueType>
void OpenMPDenseUtils::setCSRValues(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<DenseValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", size is " << numRows << " x " << numColumns )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < numRows; i++ )
    {
        // Initialize complete row with zero values

        for( IndexType j = 0; j < numColumns; ++j )
        {
            DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];

            elem = static_cast<DenseValueType>(0.0);
        }

        // fill up positions for which non-zero values are given

        for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            const IndexType j = csrJA[jj];

            DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];

            elem = static_cast<DenseValueType>( csrValues[jj] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType1,typename DenseValueType2>
void OpenMPDenseUtils::copyDenseValues(
    DenseValueType1 newValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType2 oldValues[] )
{
    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < numRows; ++i )
    {
        for( IndexType j = 0; j < numColumns; ++j )
        {
            DenseValueType1& newElem = newValues[denseindex( i, j, numRows, numColumns )];
            const DenseValueType2& oldElem = oldValues[denseindex( i, j, numRows, numColumns )];

            newElem = static_cast<DenseValueType1>( oldElem );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename RowValueType, typename DenseValueType>
void OpenMPDenseUtils::getRow(
    RowValueType rowValues[],
    const DenseValueType denseValues[],
    const IndexType irow,
    const IndexType numRows,
    const IndexType numColumns )
{
    SCAI_ASSERT_LT( irow, numRows, "illegal row index" )

    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        rowValues[j] = static_cast<RowValueType>( denseValues[irow * numColumns + j] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename DiagonalValueType,typename DenseValueType>
void OpenMPDenseUtils::getDiagonal(
    DiagonalValueType diagonalValues[],
    const IndexType numDiagonalValues,
    const DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns )
{
    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numDiagonalValues; ++i )
    {
        const DenseValueType& elem = denseValues[denseindex( i, i, numRows, numColumns )];
        diagonalValues[i] = static_cast<DiagonalValueType>( elem );
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType,typename DiagonalValueType>
void OpenMPDenseUtils::setDiagonal(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DiagonalValueType diagonalValues[],
    const IndexType numDiagonalValues )
{
    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numDiagonalValues; ++i )
    {
        DenseValueType& elem = denseValues[denseindex( i, i, numRows, numColumns )];
        elem = static_cast<DenseValueType>( diagonalValues[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
void OpenMPDenseUtils::setDiagonalValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType diagonalValue )
{
    IndexType numDiagonalValues = std::min( numRows, numColumns );

    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numDiagonalValues; ++i )
    {
        DenseValueType& elem = denseValues[denseindex( i, i, numRows, numColumns )];
        elem = diagonalValue;
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
void OpenMPDenseUtils::setValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType val )
{
    // Parallel initialization very important for efficient  allocation

    #pragma omp parallel for schedule ( SCAI_OMP_SCHEDULE )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseValues[denseindex( i, j, numRows, numColumns )] = val;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
void OpenMPDenseUtils::scaleValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType val )
{
    if ( val == scai::common::constants::ZERO )
    {
        // use setValue, as scaleValue might not work correctly on uninitialized data

        setValue( denseValues, numRows, numColumns, val );
    }
    else
    {
        #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType j = 0; j < numColumns; ++j )
            {
                DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
                elem *= val;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType, typename OtherValueType>
void OpenMPDenseUtils::scaleRows(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const OtherValueType rowValues[] )
{
    #pragma omp parallel for schedule ( SCAI_OMP_SCHEDULE )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const DenseValueType scaleValue = static_cast<DenseValueType>( rowValues[i] );

        // scale the whole row with this value

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseValues[denseindex( i, j, numRows, numColumns )] *= scaleValue;
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPDenseUtils::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // lower priority

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define KREGISTRY_DENSE2_REGISTER(z, J, TYPE )                                                                              \
    KernelRegistry::set<DenseKernelTrait::setCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setCSRValues, Host, flag );       \
    KernelRegistry::set<DenseKernelTrait::getCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getCSRValues, Host, flag );       \
    KernelRegistry::set<DenseKernelTrait::copyDenseValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( copyDenseValues, Host, flag ); \
    KernelRegistry::set<DenseKernelTrait::getDiagonal<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getDiagonal, Host, flag );         \
    KernelRegistry::set<DenseKernelTrait::setDiagonal<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setDiagonal, Host, flag );         \
    KernelRegistry::set<DenseKernelTrait::getRow<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getRow, Host, flag );                   \
    KernelRegistry::set<DenseKernelTrait::scaleRows<TYPE, ARITHMETIC_HOST_TYPE_##J> >( scaleRows, Host, flag );             \

#define KREGISTRY_DENSE_REGISTER(z, I, _)                                                                                   \
    KernelRegistry::set<DenseKernelTrait::nonZeroValues<ARITHMETIC_HOST_TYPE_##I> >( nonZeroValues, Host, flag );           \
    KernelRegistry::set<DenseKernelTrait::getCSRSizes<ARITHMETIC_HOST_TYPE_##I> >( getCSRSizes, Host, flag );               \
    KernelRegistry::set<DenseKernelTrait::setValue<ARITHMETIC_HOST_TYPE_##I> >( setValue, Host, flag );                     \
    KernelRegistry::set<DenseKernelTrait::scaleValue<ARITHMETIC_HOST_TYPE_##I> >( scaleValue, Host, flag );                 \
    KernelRegistry::set<DenseKernelTrait::setDiagonalValue<ARITHMETIC_HOST_TYPE_##I> >( setDiagonalValue, Host, flag );     \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, KREGISTRY_DENSE2_REGISTER, ARITHMETIC_HOST_TYPE_##I )                        \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, KREGISTRY_DENSE_REGISTER, _ )

#undef KREGISTRY_DENSE_REGISTER
#undef KREGISTRY_DENSE2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils::OpenMPDenseUtils()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

OpenMPDenseUtils::~OpenMPDenseUtils()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils OpenMPDenseUtils::guard;

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
