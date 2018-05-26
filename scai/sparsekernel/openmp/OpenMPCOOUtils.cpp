/**
 * @file OpenMPCOOUtils.cpp
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
 * @brief Implementation of COO storage utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/openmp/OpenMPCOOUtils.hpp>

// local library
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tracing.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/OpenMP.hpp>

#include <functional>

namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPCOOUtils::logger, "OpenMP.COOUtils" )

/* --------------------------------------------------------------------------- */

IndexType OpenMPCOOUtils::getRowStartPos( const IndexType i, 
                                          const IndexType cooIA[], 
                                          const IndexType numValues )
{
    IndexType first = 0;
    IndexType last  = numValues;

    while ( first < last )
    {
        IndexType middle = first + ( last - first ) / 2;

        if ( cooIA[middle] == i )
        {
            while ( middle > 0 && cooIA[middle - 1] == i )
            {
                middle--;
            }
            return middle;
        }
        else if ( cooIA[middle] > i )
        {
            last = middle;
        }
        else 
        {
            first = middle + 1;
        }
    }

    return invalidIndex;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCOOUtils::getValuePos( const IndexType i, const IndexType j,
                                       const IndexType cooIA[], const IndexType cooJA[],
                                       const IndexType numValues )
{
    IndexType first = 0;
    IndexType last  = numValues;

    while ( first < last )
    {
        IndexType middle = first + ( last - first ) / 2;

        if ( cooIA[middle] == i && cooJA[middle] == j )
        {
            return middle;
        }
        else if ( cooIA[middle] > i )
        {
            last = middle;
        }
        else if ( cooIA[middle] < i )
        {
            first = middle + 1;
        }
        else if ( cooJA[middle] > j )
        {
            last = middle;
        }
        else
        {
            first = middle + 1;
        }
    }

    return invalidIndex;
}

/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::offsets2ia(
    IndexType cooIA[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "build cooIA( " << numValues << " ) from csrIA( " << ( numRows + 1 ) << " )" )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        // fill all entries for row i

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            cooIA[jj] = i;
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::ia2offsets(
    IndexType csrIA[],
    const IndexType numRows,
    const IndexType cooIA[],
    const IndexType numValues )
{
    csrIA[0] = 0;

    // Be careful, this loop is inherent sequential 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offs = csrIA[i];

        while ( offs < numValues && cooIA[offs] == i )
        {
            offs++;
        }

        csrIA[i + 1] = offs;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::scaleRows(
    ValueType cooValues[],
    const ValueType rowValues[],
    const IndexType cooIA[],
    const IndexType numValues )
{
    SCAI_LOG_INFO( logger, "scaleRows in COO format" )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numValues; ++i )
    {
        cooValues[i] *= rowValues[cooIA[i]];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const common::MatrixOp op )
{
    bool sortOptimized = true;     // should only be set if COO data is sorted by rows

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger,
                       "normalGEMV<" << TypeTraits<ValueType>::id() << "> launch it asynchronously" )

        syncToken->run( std::bind( normalGEMV<ValueType>,
                                   result, alpha, x, beta, y,
                                   numRows, numColumns, numValues, 
                                   cooIA, cooJA, cooValues, op ) );
        return;
    }

    const IndexType nResult = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">,"
                   << " result[" << nResult << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    utilskernel::OpenMPUtils::binaryOpScalar( result, y, beta, nResult, common::BinaryOp::MULT, false );

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.COO.GEMV_t" )

            #pragma omp for 

            for ( IndexType k = 0; k < numValues; ++k )
            {
                IndexType i = cooIA[k];
                IndexType j = cooJA[k];
                // we must use atomic updates as different threads might update same row i
                const ValueType resultUpdate = alpha * cooValues[k] * x[i];
                // thread-safe atomic update
                atomicAdd( result[j], resultUpdate );
            }
        }
    }
    else if ( op == common::MatrixOp::NORMAL && sortOptimized )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.COO.GEMV_n" )

            // determine for each thread a range of IA so that there is no row overlap
            // by this way we can avoid the atomic add 

            IndexType k_lb;
            IndexType k_ub;

            omp_get_my_range( k_lb, k_ub, numValues );

            //      0   0   0   0   1   1   1   1  2  2   2  3  3   3
            //                          | lb                    | ub

            if ( k_lb > 0 )
            {
                // move lb up where a new row starts

                while ( k_lb < numValues && cooIA[k_lb] == cooIA[k_lb - 1] )
                {
                    k_lb++;
                }
            }

            if ( k_lb < k_ub )
            {
                // move ub up where a new row starts

                while ( k_ub < numValues && cooIA[k_ub] == cooIA[k_ub - 1] )
                {
                    k_ub++;
                }
            }

            for ( IndexType k = k_lb; k < k_ub; ++k )
            {
                result[cooIA[k]] += alpha * cooValues[k] * x[cooJA[k]];
            }
        }
    }
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.COO.GEMV_n" )
            #pragma omp for
    
            for ( IndexType k = 0; k < numValues; ++k )
            {
                IndexType i = cooIA[k];
                IndexType j = cooJA[k];
                // we must use atomic updates as different threads might update same row i
                const ValueType resultUpdate = alpha * cooValues[k] * x[j];
                atomicAdd( result[i], resultUpdate );
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "matrix op " << op << " unsupported here." )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::jacobi(
    ValueType* solution,
    const IndexType cooNumValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // solution = omega * ( rhs - B * oldSolution ) * dinv + ( 1 - omega * oldSolution
    // done in two steps
    // solution = omega * rhs * dinv + ( 1 - omega * oldSolution
    // solution -= omega * B * oldSolution * dinv

    ValueType omega1 = ValueType( 1 ) - omega;

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        ValueType diag = 0;
        ValueType vSum = rhs[i];

        IndexType k = getRowStartPos( i, cooIA, cooNumValues );

        SCAI_ASSERT_NE_ERROR( k, invalidIndex, "no entry in row " << i )

        while ( k < cooNumValues && cooIA[k] == i )
        {
            IndexType j = cooJA[k];
 
            if ( j != i )
            {
                vSum -= cooValues[k] * oldSolution[j];
            }
            else
            {
                diag = cooValues[k];
            }
            k++;
        }

        SCAI_ASSERT_NE_ERROR( diag, ValueType( 0 ), "no diagonal element or zero in row " << i )

        solution[i] = omega * vSum / diag + omega1 * oldSolution[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::getDiagonal(
    ValueType diagonal[],
    const IndexType numDiagonals,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const IndexType numValues )
{
    #pragma omp parallel for 
    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        diagonal[i] = ValueType( 0 );
    }

    #pragma omp parallel for 
    for ( IndexType k = 0; k < numValues; ++k )
    {
        IndexType i = cooIA[k];

        if ( cooJA[k] == i )
        {
            SCAI_ASSERT_VALID_INDEX_DEBUG( i, numDiagonals, "number of diagonals illegally specified" )
            diagonal[i] = cooValues[k];
        }
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPCOOUtils::hasDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const IndexType numValues )
{
    bool diagonalProperty = true;

    #pragma omp parallel for 
    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        if ( !diagonalProperty )
        {
            continue;   // no need for further checks
        }

        IndexType k = getRowStartPos( i, cooIA, numValues );

        if ( k == invalidIndex )
        {
            // row i has no entries at all
            diagonalProperty = false;
            continue;
        }

        bool found = false;

        while ( k < numValues && cooIA[k] == i )
        {
            // traverse all entrys of row

            if ( cooJA[k]  == i )
            {
                found = true;
                break;
            }
            k++;
        }
 
        if ( !found ) 
        {
            diagonalProperty = false;
        }
    }

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::setDiagonalV(
    ValueType cooValues[],
    const ValueType diagonal[],
    const IndexType numDiagonals,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const IndexType numValues )
{   
    #pragma omp parallel for 
    for ( IndexType k = 0; k < numValues; ++k )
    {   
        IndexType i = cooIA[k];
        
        if ( cooJA[k] == i )
        {   
            SCAI_ASSERT_VALID_INDEX_DEBUG( i, numDiagonals, "number of diagonals illegally specified" )
            cooValues[k] = diagonal[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::setDiagonal(
    ValueType cooValues[],
    const ValueType diagonalValue,
    const IndexType numDiagonals,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const IndexType numValues )
{   
    #pragma omp parallel for 
    for ( IndexType k = 0; k < numValues; ++k )
    {
        IndexType i = cooIA[k];

        if ( cooJA[k] == i )
        {
            SCAI_ASSERT_VALID_INDEX_DEBUG( i, numDiagonals, "number of diagonals illegally specified" )
            cooValues[k] = diagonalValue;
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCOOUtils::getRow(
    IndexType& offset,
    const IndexType cooIA[],
    const IndexType numValues,
    const IndexType i )
{
    offset = getRowStartPos( i, cooIA, numValues );

    IndexType count = 0;

    if ( offset != invalidIndex )
    {
        while ( ( offset + count ) < numValues && cooIA[offset + count] == i )
        {
            count++;
        }
    }

    return count;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCOOUtils::getColumn(
    IndexType positions[],
    const IndexType cooJA[],
    const IndexType numValues,
    const IndexType j )
{
    IndexType count = 0;

    if ( positions == NULL )
    {
        // only count the number of non-zero entries for column j

        #pragma omp parallel for reduction( +:count )
        for ( IndexType k = 0; k < numValues; ++k )
        {
            if ( cooJA[k] == j )
            {
                count++;
            }
        }
    }
    else
    { 
        // store the positions of the non-zero entries for column j

        // be careful about parallelization as we need sorted positions

        for ( IndexType k = 0; k < numValues; ++k )
        {
            if ( cooJA[k] == j )
            {
                positions[count++] = k;
            }
        }
    }

    return count;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<COOKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<COOKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<COOKernelTrait::offsets2ia>( offsets2ia, ctx, flag );
    KernelRegistry::set<COOKernelTrait::getColumn>( getColumn, ctx, flag );
    KernelRegistry::set<COOKernelTrait::getRow>( getRow, ctx, flag );
}

template<typename ValueType>
void OpenMPCOOUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<COOKernelTrait::getDiagonal<ValueType> >( getDiagonal, ctx, flag );
    KernelRegistry::set<COOKernelTrait::setDiagonal<ValueType> >( setDiagonal, ctx, flag );
    KernelRegistry::set<COOKernelTrait::setDiagonalV<ValueType> >( setDiagonalV, ctx, flag );
    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<COOKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<COOKernelTrait::scaleRows<ValueType> >( scaleRows, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils::OpenMPCOOUtils()
{
    SCAI_LOG_INFO( logger, "register COOUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPCOOUtils::~OpenMPCOOUtils()
{
    SCAI_LOG_INFO( logger, "unregister COOUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils OpenMPCOOUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
