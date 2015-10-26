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
#include <scai/lama/UtilKernelTrait.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/OpenMP.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

namespace lama
{

using std::abs;
// so we can use abs for float and double and abs for Complex<ValueType>

using common::getScalarType;

SCAI_LOG_DEF_LOGGER( OpenMPDenseUtils::logger, "OpenMP.DenseUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

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

            if( abs( value ) > eps )
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
                   "get CSRValues<" << getScalarType<DenseValueType>() << ", " << getScalarType<CSRValueType>() << ">" << ", size is " << numRows << " x " << numColumns )

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

            if( abs( value ) > eps )
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
                   "set CSRValues<" << getScalarType<DenseValueType>() << ", " << getScalarType<CSRValueType>() << ">" << ", size is " << numRows << " x " << numColumns )

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
void OpenMPDenseUtils::scaleValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType val )
{
    if( val == scai::common::constants::ZERO )
    {
        // this solution can also deal with undefined data

        #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)
        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType j = 0; j < numColumns; ++j )
            {
                DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
                elem = static_cast<DenseValueType>(0.0);
            }
        }
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
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPDenseUtils::registerKernelFunctions()
{
    using namespace scai::kregistry;

    // ctx will contain the context for which registration is done, here Host

    common::ContextType ctx = common::context::Host;

#define KREGISTRY_DENSE2_REGISTER(z, J, TYPE )                                                                          \
    /* Conversions  */                                                                                                  \
    KernelRegistry::set<DenseKernelTrait::setCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setCSRValues, ctx );       \
    KernelRegistry::set<DenseKernelTrait::getCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getCSRValues, ctx );       \
    /* Copy  */                                                                                                         \
    KernelRegistry::set<DenseKernelTrait::copyDenseValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( copyDenseValues, ctx ); \
    KernelRegistry::set<DenseKernelTrait::getDiagonal<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getDiagonal, ctx );         \
    KernelRegistry::set<DenseKernelTrait::setDiagonal<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setDiagonal, ctx );         \

#define KREGISTRY_DENSE_REGISTER(z, I, _)                                                                               \
    /* Counting  */                                                                                                     \
    KernelRegistry::set<DenseKernelTrait::getCSRSizes<ARITHMETIC_HOST_TYPE_##I> >( getCSRSizes, ctx );               \
    /* Modify  */                                                                                                       \
    KernelRegistry::set<DenseKernelTrait::setDiagonalValue<ARITHMETIC_HOST_TYPE_##I> >( setDiagonalValue, ctx );     \
    KernelRegistry::set<DenseKernelTrait::scaleValue<ARITHMETIC_HOST_TYPE_##I> >( scaleValue, ctx );                 \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, KREGISTRY_DENSE2_REGISTER, ARITHMETIC_HOST_TYPE_##I )                         \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, KREGISTRY_DENSE_REGISTER, _ )

#undef KREGISTRY_DENSE_REGISTER
#undef KREGISTRY_DENSE2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the DenseUtils routines                           */
/* --------------------------------------------------------------------------- */

bool OpenMPDenseUtils::registerInterface()
{
    registerKernelFunctions();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPDenseUtils::initialized = registerInterface();

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
