/**
 * @file OpenMPDIAUtils.cpp
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
 * @brief Implementation of DIA utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// for dll_import
#include <scai/lama/openmp/OpenMPUtils.hpp>

// local library
#include <scai/lama/UtilKernelTrait.hpp>
#include <scai/lama/openmp/OpenMPDIAUtils.hpp>
#include <scai/lama/openmp/OpenMP.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Assert.hpp>
#include <scai/common/ScalarType.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <cmath>

namespace scai
{

namespace lama
{

using std::abs;
using common::getScalarType;
using tasking::SyncToken;

SCAI_LOG_DEF_LOGGER( OpenMPDIAUtils::logger, "OpenMP.DIAUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPDIAUtils::absMaxVal(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    ValueType maxValue = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType d = 0; d < numDiagonals; ++d )
            {
                const IndexType j = i + diaOffsets[d];

                if( ( j < 0 ) || ( j >= numColumns ) )
                {
                    continue;
                }

                const ValueType val = abs( diaValues[i + d * numRows] );

                if( val > threadVal )
                {
                    threadVal = val;
                }
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if( threadVal > maxValue )
            {
                maxValue = threadVal;
            }
        }
    }

    return maxValue;
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType,typename CSRValueType>
void OpenMPDIAUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[],
    const DIAValueType eps )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << getScalarType<DIAValueType>() << ", " << getScalarType<CSRValueType>() 
                    << ">" << ", #rows = " << numRows << ", #diagonals = " << numDiagonals 
                    << ", #non-zero values = " << csrIA[numRows] << ", diagonalFlag = " << diagonalFlag )

    // we cannot check for correct sizes, but at least for valid pointers

    if( numDiagonals == 0 )
    {
        if( diagonalFlag )
        {
            IndexType n = std::min( numRows, numColumns );

            SCAI_ASSERT_EQUAL_DEBUG( n, csrIA[numRows] )

            #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType i = 0; i < n; i++ )
            {
                csrJA[i] = i;
                csrValues[i] = static_cast<CSRValueType>(0.0);
            }
        }
        else
        {
            SCAI_ASSERT_EQUAL_DEBUG( 0, csrIA[numRows] )
        }

        return;
    }

    if( numDiagonals > 0 )
    {
        SCAI_ASSERT_DEBUG( diaOffsets != NULL, "offset array of DIA data is NULL" )

        if( numRows > 0 )
        {
            SCAI_ASSERT_DEBUG( diaValues != NULL, "value array of DIA data is NULL" )
        }
    }

    // go through the DIA the same way again and copy the non-zeros

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA->CSR_values" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = csrIA[i];

            IndexType ii0 = 0; // first index of diagonal

            if( diagonalFlag && ( i < numColumns ) )
            {
                // store main diagonal at first, must be first diagonal

                SCAI_ASSERT_EQUAL_ERROR( diaOffsets[0], 0 )

                csrJA[offset] = i;
                csrValues[offset] = static_cast<CSRValueType>( diaValues[i] );

                SCAI_LOG_TRACE( logger,
                                "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                offset++;
                ii0 = 1;
            }

            for( IndexType ii = ii0; ii < numDiagonals; ii++ )
            {
                IndexType j = i + diaOffsets[ii];

                if( j < 0 )
                {
                    continue;
                }

                if( j >= numColumns )
                {
                    break;
                }

                const DIAValueType value = diaValues[i + ii * numRows];

                bool nonZero = abs( value ) > eps;

                if( nonZero )
                {
                    csrJA[offset] = j;
                    csrValues[offset] = static_cast<CSRValueType>( value );
                    SCAI_LOG_TRACE( logger,
                                    "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                    offset++;
                }
            }

            SCAI_ASSERT_EQUAL_DEBUG( offset, csrIA[i + 1] )
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType>
void OpenMPDIAUtils::getCSRSizes(
    IndexType csrSizes[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[],
    const DIAValueType eps )
{
    SCAI_LOG_INFO( logger,
                   "get CSRSizes<" << getScalarType<DIAValueType>() << "> for DIA matrix " << numRows << " x " << numColumns 
                    << ", #diagonals = " << numDiagonals << ", eps = " << eps << ", diagonalFlag = " << diagonalFlag )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; i++ )
    {
        IndexType count = 0;

        if( diagonalFlag && ( i < numColumns ) )
        {
            count = 1;
        }

        for( IndexType ii = 0; ii < numDiagonals; ii++ )
        {
            IndexType j = i + diaOffsets[ii]; // column index

            if( j < 0 )
            {
                continue;
            }

            if( j >= numColumns )
            {
                break;
            }

            bool nonZero = abs( diaValues[i + ii * numRows] ) > eps;

            if( diagonalFlag && ( i == j ) )
            {
                nonZero = false; // already counted
            }

            if( nonZero )
            {
                count++;
            }
        }

        csrSizes[i] = count;

        SCAI_LOG_TRACE( logger, "#entries in row " << i << ": " << count )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() << ">, result[" << numRows << "] = " << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << getScalarType<ValueType>() << ">, n = " << numRows << ", d = " << numDiagonals )

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    OpenMPUtils::setScale( result, beta, y, numRows );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.normalGEMV" )

        #pragma omp for schedule ( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; i++ )
        {
            ValueType accu = static_cast<ValueType>(0.0);

            for( IndexType ii = 0; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + diaOffsets[ii];

                if( j >= numColumns )
                {
                    break;
                }

                if( j >= 0 )
                {
                    accu += diaValues[ii * numRows + i] * x[j];
                }
            }

            result[i] += alpha * accu;
        }

        if( SCAI_LOG_TRACE_ON( logger ) )
        {
            std::cout << "NormalGEMV: result = ";

            for( IndexType i = 0; i < numRows; ++i )
            {
                std::cout << " " << result[i];
            }

            std::cout << std::endl;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[],
    SyncToken* syncToken )
{
    if( !syncToken )
    {
        normalGEMV( result, alpha, x, beta, y, numRows, numColumns, numDiagonals, diaOffsets, diaValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "no asynchronous support due to boost problem" )

        // asynchronous execution is done by calling an own thread at higher level
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[],
    SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "normalGEVM<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() << ">, result[" << numRows << "] = " << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    SCAI_LOG_INFO( logger,
                   "normalGEVM<" << getScalarType<ValueType>() << ">, n = " << numRows << ", d = " << numDiagonals )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    OpenMPUtils::setScale( result, beta, y, numColumns );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.normalGEVM" )

        #pragma omp for schedule ( SCAI_OMP_SCHEDULE )

        for( IndexType k = 0; k < numColumns; ++k )
        {
            ValueType accu = static_cast<ValueType>(0.0);

            for( IndexType i = 0; i < numRows; i++ )
            {
                for( IndexType ii = 0; ii < numDiagonals; ++ii )
                {
                    const IndexType j = i + diaOffsets[ii];

                    if( j >= numColumns )
                    {
                        break;
                    }

                    if( j == k )
                    {
                        accu += diaValues[ii * numRows + i] * x[i];
                    }
                }
            }

            result[k] += alpha * accu;
        }
    }

    if( SCAI_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEVM: result = ";

        for( IndexType i = 0; i < numRows; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

/* --------------------------------------------------------------------------- */
/*  Jacobi                                                                     */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::jacobi(
    ValueType solution[],
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffset[],
    const ValueType diaValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows,
    class SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", #cols = " << numColumns << ", #diagonals = " << numDiagonals << ", omega = " << omega )

    SCAI_ASSERT_EQUAL_DEBUG( 0, diaOffset[0] )
    // main diagonal must be first

    if( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.Jacobi" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            ValueType diag = diaValues[i]; // diagonal is first

            for( IndexType ii = 1; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + diaOffset[ii];

                if( j >= numColumns )
                {
                    break;
                }

                if( j >= 0 )
                {
                    temp -= diaValues[ii * numRows + i] * oldSolution[j];
                }
            }

            solution[i] = omega * ( temp / diag ) + (static_cast<ValueType>(1.0) - omega ) * oldSolution[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPDIAUtils::registerKernels()
{
    using kregistry::KernelRegistry;

    // ctx will contain the context for which registration is done, here Host

    common::ContextType ctx = common::context::Host;

    // use of BOOST_PP_REPEAT to register for all value types
    // use of nested BOOST_PP_REPEAT to get all pairs of value types for conversions

#define LAMA_DIA_UTILS2_REGISTER(z, J, TYPE )                                                                    \
    KernelRegistry::set<DIAKernelTrait::getCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getCSRValues, ctx );  \

#define LAMA_DIA_UTILS_REGISTER(z, I, _)                                                                 \
    KernelRegistry::set<DIAKernelTrait::getCSRSizes<ARITHMETIC_HOST_TYPE_##I> >( getCSRSizes, ctx );  \
    KernelRegistry::set<DIAKernelTrait::absMaxVal<ARITHMETIC_HOST_TYPE_##I> >( absMaxVal, ctx );      \
    KernelRegistry::set<DIAKernelTrait::normalGEMV<ARITHMETIC_HOST_TYPE_##I> >( normalGEMV, ctx );    \
    KernelRegistry::set<DIAKernelTrait::normalGEVM<ARITHMETIC_HOST_TYPE_##I> >( normalGEVM, ctx );    \
    KernelRegistry::set<DIAKernelTrait::jacobi<ARITHMETIC_HOST_TYPE_##I> >( jacobi, ctx );            \
                                                                                                         \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT,                                                           \
                     LAMA_DIA_UTILS2_REGISTER,                                                           \
                     ARITHMETIC_HOST_TYPE_##I )

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_DIA_UTILS_REGISTER, _ )

#undef LAMA_DIA_UTILS_REGISTER
#undef LAMA_DIA_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the DIAUtils routines                             */
/* --------------------------------------------------------------------------- */

bool OpenMPDIAUtils::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPDIAUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
