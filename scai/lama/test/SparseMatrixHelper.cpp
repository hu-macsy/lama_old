/**
 * @file SparseMatrixHelper.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief SparseMatrixHelper.cpp
 * @author Lauretta Schubert
 * @date 20.04.2012
 */

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama/test/SparseMatrixHelper.hpp>

template<typename ValueType>
bool compareMatrices( const scai::lama::CSRSparseMatrix<ValueType>& matrix, const scai::lama::CSRSparseMatrix<ValueType>& other )
{
    IndexType nnu = matrix.getLocalNumRows();

    if ( !( nnu == other.getLocalNumRows() && matrix.getLocalNumColumns() == other.getLocalNumColumns()
            && matrix.getLocalNumValues() == other.getLocalNumValues() ) )
    {
        return false;
    }

    // Note: getLocalData returns CSRStorage<ValueType>, so we can use getIA
    const scai::lama::CSRStorage<ValueType>& matrixLocal = matrix.getLocalStorage();
    const scai::lama::CSRStorage<ValueType>& otherLocal = other.getLocalStorage();
    scai::hmemo::ReadAccess<IndexType> ia( matrixLocal.getIA() );
    scai::hmemo::ReadAccess<IndexType> iaOther( otherLocal.getIA() );

    for ( IndexType i = 0; i < nnu + 1; ++i )
    {
        if ( ia[i] != iaOther[i] )
        {
            return false;
        }
    }

    scai::hmemo::ReadAccess<IndexType> ja( matrixLocal.getJA() );
    scai::hmemo::ReadAccess<IndexType> jaOther( otherLocal.getJA() );
    scai::hmemo::ReadAccess<ValueType> data( matrixLocal.getValues() );
    scai::hmemo::ReadAccess<ValueType> dataOther( otherLocal.getValues() );

    for ( IndexType i = 0; i < nnu; ++i )
    {
        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            IndexType j = ja[jj];
            bool found = false;

            for ( IndexType kk = iaOther[i]; kk < iaOther[i + 1]; ++kk )
            {
                if ( j == jaOther[kk] )
                {
                    found = true;

                    if ( data[jj] != dataOther[kk] )
                    {
                        return false;
                    }
                }
            }

            if ( !found )
            {
                return false;
            }
        }
    }

    return true;
}

#define COMPARE_INSTANTIATE( type )   \
    template  \
    COMMON_DLL_IMPORTEXPORT \
    bool compareMatrices( const scai::lama::CSRSparseMatrix<type>& matrix, const scai::lama::CSRSparseMatrix<type>& other );

SCAI_COMMON_LOOP( COMPARE_INSTANTIATE, SCAI_ARITHMETIC_HOST )
