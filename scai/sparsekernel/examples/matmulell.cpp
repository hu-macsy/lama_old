/**
 * @file sparsekernel/examples/matmulell.cpp
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
 * @brief Matrix vector multiplication using the ELL-format
 * @author Eric Schricker
 * @date 21.03.2016
 */

#include <scai/common/ContextType.hpp>

#include <scai/sparsekernel/ELLKernelTrait.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace common;
using namespace kregistry;
using namespace sparsekernel;

/*
 * initialize array with a given value
 */
template<typename ValueType>
static void init( const IndexType n, ValueType* x, const ValueType val )
{
    for ( IndexType i = 0; i < n; ++i )
    {
        x[i] = val;
    }
}

/*
 * print array
 */
template<typename ValueType>
static void print( const IndexType n, const ValueType* x )
{
    for ( IndexType i = 0; i < n; ++i )
    {
        std::cout << x[i] << " ";
    }

    std::cout << std::endl;
}

template<typename ValueType>
static void multiplication()
{
    IndexType m, n, max_nnz;
    m = 4;
    n = 3;
    max_nnz = 2;
    const ValueType one = ValueType( 1 );
    const ValueType zero = ValueType( 0 );
    /* Matrix
     *
     * Dense:
     *
     *   7  |  1.2  | 0
     *   0  |   3   | 0
     *   4  |   0   | 1
     *   0  |  7.5  | 0
     *
     *  ELL:
     *
     *       values:               ja:          sizes:
     *   7  |  1.2  | 0           0 | 1           2
     *   3  |   0   | 0           1 | 0           1
     *   4  |   1   | 0           0 | 2           2
     *  7.5 |   0   | 0           1 | 0           1
     */
    // allocate memory

    ValueType* values = new ValueType[ m * max_nnz ];
    IndexType* ja = new IndexType[ m * max_nnz ];
    IndexType* sizes = new IndexType[ m ];
    // init
    init( m * max_nnz, values, zero );
    init( m * max_nnz, ja, IndexType( 0 ) );
    init( m, sizes, IndexType( 0 ) );
    // Values-Array
    //     i             j  Index
    values[0 + m * 0] = 7;
    values[1 + m * 0] = 3;
    values[2 + m * 0] = 4;
    values[3 + m * 0] = 7.5;
    values[0 + m * 1] = 1.2;
    values[2 + m * 1] = 1;
    // JA-Array
    // i             j
    ja[0 + m * 0] = 0;
    ja[1 + m * 0] = 1;
    ja[2 + m * 0] = 0;
    ja[3 + m * 0] = 1;
    ja[0 + m * 1] = 1;
    ja[2 + m * 1] = 2;
    // Sizes-Array
    //    i
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 2;
    sizes[3] = 1;
    // x-Vector
    ValueType* x = new ValueType[ n ];
    init( n, x, one );
    // b-Vector
    ValueType* b = new ValueType[ m ];
    init( m, b, zero );
    // Get Function from Registry
    KernelTraitContextFunction<ELLKernelTrait::normalGEMV<ValueType> > gemv;
    gemv[ContextType::Host]( b, one, x, zero, b, m, n, max_nnz, sizes, ja, values, common::MatrixOp::NORMAL );
    std::cout << "Vector b: ";
    print( m, b );
    // free memory
    delete[] values;
    delete[] ja;
    delete[] sizes;
    delete[] x;
    delete[] b;
}

int main()
{
    std::cout << "float --> ";
    multiplication<float>();
    std::cout << "double --> ";
    multiplication<double>();
    std::cout << "long double --> ";
    multiplication<long double>();
#ifdef SCAI_COMPLEX_SUPPORTED
    std::cout << "ComplexFloat --> ";
    multiplication<ComplexFloat>();
    std::cout << "ComplexDouble --> ";
    multiplication<ComplexDouble>();
    std::cout << "ComplexLongDouble --> ";
    multiplication<ComplexLongDouble>();
#endif
}


