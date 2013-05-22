/**
 * @file SparseMatrixHelper.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief SparseMatrixHelper.cpp
 * @author Lauretta Schubert
 * @date 20.04.2012
 * $Id$
 */

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/HostReadAccess.hpp>

#include <test/SparseMatrixHelper.hpp>

template<typename ValueType>
bool compareMatrices( const lama::CSRSparseMatrix<ValueType>& matrix, const lama::CSRSparseMatrix<ValueType>& other )
{
    lama::IndexType nnu = matrix.getLocalNumRows();

    if ( !( nnu == other.getLocalNumRows() && matrix.getLocalNumColumns() == other.getLocalNumColumns()
            && matrix.getLocalNumValues() == other.getLocalNumValues() ) )
    {
        return false;
    }

    // Note: getLocalData returns CSRStorage<ValueType>, so we can use getIA

    const lama::CSRStorage<ValueType>& matrixLocal = matrix.getLocalStorage();
    const lama::CSRStorage<ValueType>& otherLocal = other.getLocalStorage();

    lama::HostReadAccess<lama::IndexType> ia( matrixLocal.getIA() );
    lama::HostReadAccess<lama::IndexType> iaOther( otherLocal.getIA() );

    for ( lama::IndexType i = 0; i < nnu + 1; ++i )
    {
        if ( ia[i] != iaOther[i] )
        {
            return false;
        }
    }

    lama::HostReadAccess<lama::IndexType> ja( matrixLocal.getJA() );
    lama::HostReadAccess<lama::IndexType> jaOther( otherLocal.getJA() );
    lama::HostReadAccess<ValueType> data( matrixLocal.getValues() );
    lama::HostReadAccess<ValueType> dataOther( otherLocal.getValues() );

    for ( lama::IndexType i = 0; i < nnu; ++i )
    {
        for ( lama::IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            lama::IndexType j = ja[jj];
            bool found = false;
            for ( lama::IndexType kk = iaOther[i]; kk < iaOther[i + 1]; ++kk )
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

template
#ifdef WIN32
__declspec(dllexport)
#endif
bool compareMatrices( const lama::CSRSparseMatrix<double> & matrix, const lama::CSRSparseMatrix<double> & other );

template
#ifdef WIN32
__declspec(dllexport)
#endif
bool compareMatrices( const lama::CSRSparseMatrix<float> & matrix, const lama::CSRSparseMatrix<float> & other );

