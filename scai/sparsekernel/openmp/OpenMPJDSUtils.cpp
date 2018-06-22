/**
 * @file OpenMPJDSUtils.cpp
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
 * @brief Implementation of JDS utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/openmp/OpenMPJDSUtils.hpp>

// local library
#include <scai/sparsekernel/JDSKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tracing.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <memory>
#include <functional>

namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( OpenMPJDSUtils::logger, "OpenMP.JDSUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPJDSUtils::getRow(
    ValueType row[],
    const IndexType i,
    const IndexType numColumns,
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType ja[],
    const ValueType values[] )
{
    SCAI_REGION( "OpenMP.JDS.getRow" )

    SCAI_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        row[j] = ValueType( 0 );
    }

    IndexType ii;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    IndexType k = 0;

    for ( IndexType jj = 0; jj < ilg[ii]; ++jj )
    {
        row[ja[ii + k]] = values[ii + k];
        k += dlg[jj];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPJDSUtils::setRow(
    ValueType values[],
    const IndexType i,
    const IndexType numColumns,
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType ja[],
    const ValueType row[],
    const common::BinaryOp op )
{
    SCAI_REGION( "OpenMP.JDS.setRow" )

    SCAI_LOG_INFO( logger, "setRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

    IndexType ii;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    SCAI_ASSERT_VALID_INDEX_ERROR( ii, numRows, "row " << i << " not in permutuation" )

    IndexType k = ii;

    for ( IndexType jj = 0; jj < ilg[ii]; ++jj )
    {
        values[k] = applyBinary( values[k], op, static_cast<ValueType>( row[ja[k]] ) );
        k += dlg[jj];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPJDSUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType perm[],
    const IndexType ja[] )
{
    IndexType ii;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    if ( ii == numRows )
    {
        COMMON_THROWEXCEPTION( "row index " << i << " not found in perm array" )
    }

    SCAI_LOG_TRACE( logger, "row " << i << " is now " << ii << ", has " << ilg[ii] << " elements" )

    // search in the found row

    IndexType k = 0;

    IndexType pos = invalidIndex;

    for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
    {
        if ( ja[ii + k] == j )
        {
            pos = ii + k;
            break;
        }

        k += dlg[jj];
    }

    return pos;
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPJDSUtils::getDiagonalPositions(
    IndexType diagonalPositions[],
    const IndexType numDiagonals,
    const IndexType numRows,
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType perm[],
    const IndexType ja[] )
{
    #pragma omp parallel for
    for ( IndexType i = 0; i < numDiagonals; i++ )
    {
        diagonalPositions[i] = invalidIndex;
    }

    IndexType numDiagonalsFound = 0;  

    #pragma omp parallel for reduction( + : numDiagonalsFound )
    for ( IndexType ii = 0; ii < numRows; ii++ )
    {
        IndexType i = perm[ii];    // original row

        IndexType k = 0;

        SCAI_LOG_TRACE( logger, "find diag entry for row " << i << ", has " << ilg[ii] << " entries" )

        for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
        {
            SCAI_LOG_TRACE( logger, "check entry " << jj << " of " << ilg[ii] << ", j = " << ja[ii+k] )

            if ( ja[ii + k] == i )
            {
                SCAI_ASSERT_LT_ERROR( i, numDiagonals, "more diagonals than expected" )
                diagonalPositions[i] = ii + k;
                numDiagonalsFound++;
                break;
            }
    
            k += dlg[jj];
        }
    }

    return numDiagonalsFound;
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPJDSUtils::getColumnPositions(
    IndexType row[],
    IndexType pos[],
    const IndexType j,
    const IndexType numRows,
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType perm[],
    const IndexType ja[] )
{
    SCAI_REGION( "OpenMP.JDSUtils.getColumnPositions" )

    IndexType cnt  = 0;   // counts number of available row entries in column j

    #pragma omp parallel for

    for ( IndexType ii = 0; ii < numRows; ++ii )
    {
        IndexType k = 0;

        for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
        {
            IndexType p = ii + k;

            if ( ja[p] == j )
            {
                IndexType n = atomicInc( cnt );
                row[n] = perm[ii];
                pos[n] = p;
                break;
            }

            k += dlg[jj];
        }
    }

    return cnt;
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPJDSUtils::getRowPositions(
    IndexType pos[],
    const IndexType i,
    const IndexType numRows,
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType perm[] )
{
    SCAI_REGION( "OpenMP.JDSUtils.getRowPositions" )

    IndexType ii = invalidIndex;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    SCAI_ASSERT_NE_ERROR( ii, numRows, "row " << i << " not found in perm array" );

    IndexType cnt = ilg[ii];

    IndexType k = 0;

    for ( IndexType jj = 0; jj < cnt; ++jj )
    {
        pos[jj] = ii + k;
        k += dlg[jj];
        SCAI_LOG_TRACE( logger, "pos[" << jj << "] = " << pos[jj] << ", dlg[" << jj << "] = " << dlg[jj] )
    }

    return cnt;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPJDSUtils::setRows(
    ValueType jdsValues[],
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const ValueType rowValues[],
    const common::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "setRows with numRows = " << numRows )

    // Due to false sharing, use of OpenMP is not recommended here

    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType offset = i;

        ValueType rowScale = rowValues[perm[i]];

        for ( IndexType jj = 0; jj < ilg[i]; jj++ )
        {
            jdsValues[offset] = common::applyBinary( jdsValues[offset], op, rowScale );
            offset += dlg[jj];
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

    if ( numDiagonals == 0 )
    {
        return 0;
    }

    SCAI_ASSERT_EQUAL_DEBUG( numDiagonals, ilg[0] )
    // Entries in dlg filled every time there is a change in values of consecutive elements of ilg
    //
    //   i:     0  1  2  3  4  5
    // ilg:     5  5  3  3  3  1
    // nd1:     5  5  3  3  3  1
    // nd2:     5  3  3  3  1  0
    //             x        x  x
    //             |        |  |->    6
    //             |        |---->       5  5
    //             |------------->             2   2
    // dlg:                           6  5  5  2   2
    IndexType numTotal = 0;
    #pragma omp parallel for reduction( +:numTotal )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType nd1 = ilg[i];
        IndexType nd2 = 0;

        if ( i + 1 < numRows )
        {
            nd2 = ilg[i + 1];
        }

        // fill in dlg only if nd2 < nd1

        for ( IndexType j = nd2; j < nd1; j++ )
        {
            dlg[j] = i + 1;
        }

        numTotal += nd1;
    }

    return numTotal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::getCSRValues(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType jdsInversePerm[],
    const IndexType jdsILG[],
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[] )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", #values = " << csrIA[numRows] )
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.JDS.getCSR" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType ii = jdsInversePerm[i]; // where to find row i in JDS storage
            const IndexType numValuesInRow = jdsILG[ii];
            IndexType jdsOffset = ii; // run through input JDS data
            IndexType offset = csrIA[i]; // run through output data

            for ( IndexType jj = 0; jj < numValuesInRow; jj++ )
            {
                csrJA[offset + jj] = jdsJA[jdsOffset];
                csrValues[offset + jj] = jdsValues[jdsOffset];
                jdsOffset += jdsDLG[jj];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::setCSRValues(
    IndexType jdsJA[],
    ValueType jdsValues[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType /* ndlg */,
    const IndexType jdsDLG[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", #values = " << csrIA[numRows] )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.JDS.setCSR" )
        #pragma omp for 

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            IndexType i = jdsPerm[ii];
            IndexType offset = ii;

            for ( IndexType jdsJJ = 0, csrJJ = csrIA[i]; jdsJJ < jdsILG[ii]; jdsJJ++, csrJJ++ )
            {
                jdsJA[offset] = csrJA[csrJJ];
                jdsValues[offset] = csrValues[csrJJ];
                offset += jdsDLG[jdsJJ]; // index for next value of the row
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType perm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const common::MatrixOp op )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( std::bind( normalGEMV<ValueType>,
                                   result,
                                   alpha, x, beta, y,
                                   numRows, numColumns,
                                   perm, jdsILG, ndlg, jdsDLG, 
                                   jdsJA, jdsValues, op ) );
        return;
    }

    const IndexType nResult = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, op = " << op << ", result[" << nResult 
                   << "] = " << alpha << " * A( jds, ndlg = " << ndlg << " ) * x + " << beta << " * y " )

    // z = alpha * JDS * x + beta * y, remains: z += alpha * JDS * x

    if ( beta != common::Constants::ONE || &result != &y )
    {
        utilskernel::OpenMPUtils::binaryOpScalar( result, y, beta, nResult, common::BinaryOp::MULT, false );
    }

    if ( ndlg == 0 )
    {
        return; // definitively empty matrix
    }

    IndexType nonEmptyRows = jdsDLG[0];  // stands exactly for number of non-empty rows

    SCAI_LOG_DEBUG( logger, "result += alpha * A * x, #non-empty row = " << nonEmptyRows )

    #pragma omp parallel
    {
        if ( common::isTranspose( op ) )
        {
            SCAI_REGION( "OpenMP.JDS.gemv_t" )

            #pragma omp for

            for ( IndexType ii = 0; ii < nonEmptyRows; ii++ )
            {
                IndexType offset = ii;

                const ValueType tmpX   = x[perm[ii]];

                for ( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
                {
                    IndexType j = jdsJA[offset];
                    ValueType v = alpha * jdsValues[offset] * tmpX;
                    atomicAdd( result[j], v );
                    offset += jdsDLG[jj];      // jump to next value for this row
                }
            }
        }
        else
        {
            SCAI_REGION( "OpenMP.JDS.gemv_n" )
            #pragma omp for 
    
            for ( IndexType ii = 0; ii < nonEmptyRows; ii++ )
            {
                ValueType value = 0; // sums up final value
                IndexType offset = ii;
    
                for ( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
                {
                    IndexType j = jdsJA[offset];

                    value += jdsValues[offset] * x[j];
                    offset += jdsDLG[jj]; // there is next value for this row
                }
    
                // scattering needs no synchronization as values of perm are unique
                result[perm[ii]] += alpha * value;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Jacobi                                                                  */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType jdsNumDiagonals,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        // run this method with exactly this arguments by an own thread

        syncToken->run( std::bind( jacobi<ValueType>,
                                   solution, numRows, 
                                   jdsPerm, jdsILG, jdsNumDiagonals, jdsDLG,
                                   jdsJA, jdsValues, oldSolution, rhs, omega ) );
        return;
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.JDS.jacobi" )
        #pragma omp for 

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            const IndexType i = jdsPerm[ii]; // original row index
            ValueType temp = rhs[i];
            IndexType pos = ii; // index for jdsValues
            ValueType diag = 0;

            for ( IndexType j = 0; j < jdsILG[ii]; j++ )
            {
                if ( jdsJA[pos] == i )
                {
                    diag = jdsValues[pos];
                }
                else
                {
                    temp -= jdsValues[pos] * oldSolution[jdsJA[pos]];
                }

                pos += jdsDLG[j];
            }

            if ( omega == scai::common::Constants::ONE )
            {
                solution[i] = temp / diag;
            }
            else if ( 0.5 == omega )
            {
                solution[i] = omega * ( temp / diag + oldSolution[i] );
            }
            else
            {
                solution[i] = omega * ( temp / diag ) + ( static_cast<ValueType>( 1.0 ) - omega ) * oldSolution[i];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::jacobiHalo(
    ValueType solution[],
    const IndexType numRows,
    const ValueType localDiagonal[],
    const IndexType numDiagonals,
    const IndexType jdsHaloPerm[],
    const IndexType jdsHaloILG[],
    const IndexType jdsHaloDLG[],
    const IndexType jdsHaloJA[],
    const ValueType jdsHaloValues[],
    const ValueType oldSolution[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobiHalo<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    if ( numRows == 0 )
    {
        return;
    }

    if ( numDiagonals == 0 )
    {
        return;
    }

    // JDS has no row indexes, but number of non-zero rows is known
    const IndexType numNonEmptyRows = jdsHaloDLG[0];
    SCAI_LOG_DEBUG( logger, "#non empty rows = " << numNonEmptyRows )
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.JDS.jacobiHalo" )
        #pragma omp for 

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );
            const IndexType i = jdsHaloPerm[ii];
            const ValueType diag = localDiagonal[i];
            IndexType pos = ii;

            for ( IndexType j = 0; j < jdsHaloILG[ii]; j++ )
            {
                temp += jdsHaloValues[pos] * oldSolution[jdsHaloJA[pos]];
                pos += jdsHaloDLG[j];
            }

            SCAI_LOG_TRACE( logger,
                            "jds row " << ii << ", is row " << i << " in halo" << ", diag = " << diag << ", temp = " << temp )
            solution[i] -= temp * omega / diag;
            SCAI_LOG_TRACE( logger, "solution[" << i << "] = " << solution[i] )
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPJDSUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register JDSUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<JDSKernelTrait::ilg2dlg>( ilg2dlg, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getColumnPositions>( getColumnPositions, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getRowPositions>( getRowPositions, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getDiagonalPositions>( getDiagonalPositions, ctx, flag );
}

template<typename ValueType>
void OpenMPJDSUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register JDSUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<JDSKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getRow<ValueType> >( getRow, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::setRows<ValueType> >( setRows, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::setRow<ValueType> >( setRow, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::setCSRValues<ValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getCSRValues<ValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPJDSUtils::OpenMPJDSUtils()
{
    SCAI_LOG_INFO( logger, "register JDSUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPJDSUtils::~OpenMPJDSUtils()
{
    SCAI_LOG_INFO( logger, "unregister JDSUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPJDSUtils OpenMPJDSUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
