/**
 * @file blaskernel/examples/CG_BLAS.cpp
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
 * @brief Conjugated Gradient implementation based on BLAS calls.
 * @author Eric Schricker
 * @date 28.01.2016
 */

#include <scai/common/ContextType.hpp>

#include <scai/hmemo.hpp>

#include <scai/kregistry.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <iostream>

using scai::IndexType;
using scai::common::ContextType;
using scai::common::MatrixOp;
using scai::hmemo::HArray;
using scai::hmemo::ReadAccess;
using scai::hmemo::WriteAccess;
using scai::hmemo::WriteOnlyAccess;
using scai::hmemo::Context;
using scai::hmemo::ContextPtr;
using scai::kregistry::KernelTraitContextFunction;
using scai::blaskernel::BLASKernelTrait;

template<typename ValueType>
void cg( const IndexType m, const IndexType n, const IndexType lda,
         HArray<ValueType>& x, HArray<ValueType>& A, const HArray<ValueType>& b,
         const ValueType tol, const IndexType max_iter, const ContextType loc )
{
    // Needed temporary vectors
    HArray<ValueType> r( n );
    HArray<ValueType> d( n );
    HArray<ValueType> z( n );
    // Needed temporary scalars
    ValueType rOld, rNew, alpha, beta, norm;
    ContextPtr ctx = Context::getContextPtr( loc );
    // Access to vectors on used device
    ReadAccess<ValueType> A_( A, ctx );
    ReadAccess<ValueType> b_( b, ctx );
    WriteAccess<ValueType> r_( r, ctx );
    WriteAccess<ValueType> d_( d, ctx );
    WriteAccess<ValueType> z_( z, ctx );
    WriteAccess<ValueType> x_( x, ctx );
    // Used blas functions
    static KernelTraitContextFunction<BLASKernelTrait::gemv<ValueType> > gemv;
    static KernelTraitContextFunction<BLASKernelTrait::axpy<ValueType> > axpy;
    static KernelTraitContextFunction<BLASKernelTrait::copy<ValueType> > copy;
    static KernelTraitContextFunction<BLASKernelTrait::nrm2<ValueType> > nrm2;
    static KernelTraitContextFunction<BLASKernelTrait::scal<ValueType> > scal;
    static KernelTraitContextFunction<BLASKernelTrait::dot<ValueType> > dot;
    // CG init
    // r = -1 * A * x
    gemv[loc]( MatrixOp::NORMAL, m, n, -1.0, A_.get(), lda, x_.get(),
               1, 0.0, r_.get(), 1 );
    // r = b + r
    axpy[loc]( m, 1.0, b_.get(), 1, r_.get(), 1 );
    // d = r
    copy[loc]( m, r_.get(), 1, d_.get(), 1 );
    // rOld = r * r
    rOld = dot[loc]( m, r_.get(), 1, r_.get(), 1 );
    // norm = ||r||2
    norm = nrm2[loc]( m, r_.get(), 1 );

    // CG iterations
    for ( IndexType k = 0; norm > tol && k < max_iter; ++k )
    {
        // z = A * d
        gemv[loc]( MatrixOp::NORMAL, m, n, 1.0, A_.get(), lda,
                   d_.get(), 1, 0.0, z_.get(), 1 );
        // alpha = rOld / d * z
        alpha = dot[loc]( m, d_.get(), 1, z_.get(), 1 );
        alpha = rOld / alpha;
        // x = alpha * d + x
        axpy[loc]( m, alpha, d_.get(), 1, x_.get(), 1 );
        // r = -alpha * z + r
        axpy[loc]( m, -alpha, z_.get(), 1, r_.get(), 1 );
        // rNew = r * r
        rNew = dot[loc]( m, r_.get(), 1, r_.get(), 1 );
        beta = rNew / rOld;
        // d = r + beta * d
        //     d = beta * d
        scal[loc]( m, beta, d_.get(), 1 );
        //     d = r + d
        axpy[loc]( m, 1.0, r_.get(), 1, d_.get(), 1 );
        rOld = rNew;
        // norm = ||r||2
        norm = nrm2[loc]( m, r_.get(), 1 );

        if ( k >= 10 )
        {
            break;
        }
    }
}

int main()
{
    typedef double ValueType;
    const IndexType m = 2;
    // Maximum number of iterations for solver
    const IndexType max_iter = 20;
    // Pointer used for initializing data structures
    ContextPtr host = Context::getHostPtr();
    // Context used for computation
    ContextType loc = ContextType::Host;
    HArray<ValueType> A( m * m );
    {
        // Init Matrix
        WriteOnlyAccess<ValueType> A_( A, host );
        A_[0 * m + 0] = -2;
        A_[0 * m + 1] = 1;
        A_[1 * m + 0] = 1;
        A_[1 * m + 1] = 12;
    }
    HArray<ValueType> b( m );
    {
        // Init right hand side
        WriteOnlyAccess<ValueType> b_( b, host );
        b_[0] = 19;
        b_[1] = 19.017;
    }
    HArray<ValueType> x( m );
    {
        // Init solution vector with zero
        WriteOnlyAccess<ValueType> x_( x, host );
        x_[0] = 0;
        x_[1] = 0;
    }
    // Call CG
    cg( m, m, m, x, A, b, 1e-7, max_iter, loc );
    {
        // Output result vector
        ReadAccess<ValueType> x_( x );

        for ( IndexType i = 0; i < m; ++i )
        {
            std::cout << x_[i] << std::endl;
        }
    }
}

