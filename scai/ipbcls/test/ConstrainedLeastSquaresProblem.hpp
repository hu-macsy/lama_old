/**
 * @file ConstrainedLeastSquaresProblem.hpp
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
 * @brief Generation of constrained least squares problem
 * @author Andreas Borgen Longva
 * @date 12.09.2017
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

template <typename Scalar>
struct ConstrainedLeastSquaresProblem
{
    scai::lama::CSRSparseMatrix<Scalar> A;
    scai::lama::DenseVector<Scalar> x;
    scai::lama::DenseVector<Scalar> b;
    scai::lama::DenseVector<Scalar> l;
    scai::lama::DenseVector<Scalar> u;
    scai::lama::DenseVector<Scalar> lambda;
    scai::lama::DenseVector<Scalar> mu;

};

namespace detail
{
    /**
     * Generate a vector of the specified size with random elements.
     * The bias factor determines the expected proportion of non-zeros relative to zeros. That is,
     * a bias of 0.0 means that all elements will be zero, 1.0 means that all elements will be non-zero.
     */
    template <typename Scalar>
    scai::lama::DenseVector<Scalar> randomBiasedVector(scai::IndexType size, double bias);
    
    /**
     * Generate a vector of the specified size with random, non-negative elements.
     * The bias factor determines the expected proportion of non-zeros relative to zeros. That is,
     * a bias of 0.0 means that all elements will be zero, 1.0 means that all elements will be non-zero.
     */
    template <typename Scalar>
    scai::lama::DenseVector<Scalar> randomBiasedNonNegativeVector(scai::IndexType size, double bias);
    
    template <typename Scalar>
    std::vector<scai::IndexType> indicesOfNonZeros(const scai::lama::DenseVector<Scalar> & v);
}

template <typename Scalar>
ConstrainedLeastSquaresProblem<Scalar> generateConstrainedLeastSquaresProblem(scai::IndexType rows, scai::IndexType cols, double density = 0.05)
{
    using namespace scai::lama;

    const auto zeroVector = fill<DenseVector<Scalar>>(cols, 0 );
    
    ConstrainedLeastSquaresProblem<Scalar> p;
    
    p.A = zero<CSRSparseMatrix<Scalar>>( rows, cols );
    MatrixCreator::fillRandom(p.A, density);
    
    const DenseVector<Scalar> r = detail::randomBiasedNonNegativeVector<Scalar>(rows, 0.5);
//     const auto r_nonzero_indices = detail::indicesOfNonZeros(r);
    
    // TODO: Modify matrix to ensure some cases where the gradient of the unconstrained problem is identically zero
    
    p.x = detail::randomBiasedVector<Scalar>(cols, 0.9);
    p.b = p.A * p.x - r;
    p.l = p.x - detail::randomBiasedNonNegativeVector<Scalar>(cols, 0.9);
    p.u = p.x + detail::randomBiasedNonNegativeVector<Scalar>(cols, 0.9);
    
    // Add epsilons to make sure l and u do not coincide in any components

    Scalar eps = 1e-8;
    p.l -= eps;
    p.u += eps;
    
    // Gradient of the unconstrained problem
    const auto g = eval<DenseVector<Scalar>>( 2 * transpose( p.A) * r );
    
    // Lagrange multipliers for the bound constraints. Lambda is the multiplier for the lower bound,
    // while mu is the multiplier for the upper bound

    Scalar zero = 0;
    p.lambda = fill<DenseVector<Scalar>>( cols, 0 );
    p.mu = fill<DenseVector<Scalar>>(cols, 0 );
    
    // We build the matrix locally, so we permit ourselves to work with global indices directly
    for (scai::IndexType i = 0; i < cols; ++i)
    {
        const Scalar x_i = p.x.getValue(i);
        const Scalar g_i = g.getValue(i);
        
        if ( g_i >= zero)
        {
            p.mu.setValue(i, zero);
            p.lambda.setValue(i, g_i);
        } else {
            p.lambda.setValue(i, zero);
            p.mu.setValue(i, - g_i);
        }
        
        if (p.lambda.getValue(i) > zero)
        {
            p.l.setValue(i, x_i);
        }
        
        if (p.mu.getValue(i) > zero)
        {
            p.u.setValue(i, x_i);
        }
    }
    
    using scai::common::CompareOp;

    // Assert KKT conditions
    const auto rel_lower_bound = eval<DenseVector<Scalar>>( p.x - p.l );
    const auto rel_upper_bound = eval<DenseVector<Scalar>>( p.u - p.x );

    DenseVector<Scalar> lagrangian_gradient ( g );
    lagrangian_gradient -= p.lambda;
    lagrangian_gradient += p.mu;
    
    SCAI_ASSERT_ERROR(p.lambda.all( CompareOp::GE, zeroVector ), "Lagrange parameter lambda should be non-negative.");
    SCAI_ASSERT_ERROR(p.mu.all( CompareOp::GE, zeroVector ), "Lagrange parameter mu should be non-negative.");
    SCAI_ASSERT_ERROR(p.x.all( CompareOp::GE, p.l ), "x is not within lower bound constraints.");
    SCAI_ASSERT_ERROR(p.x.all( CompareOp::LE, p.u ), "x is not within upper bound constraints.");   
    SCAI_ASSERT_ERROR(lagrangian_gradient.all(CompareOp::LE, zeroVector)
                   && lagrangian_gradient.all(CompareOp::GE, zeroVector), "Lagrangian gradient is not zero");
    
    // Due to the lack of an equality operator for vectors, we must once again use both LE and GE

    SCAI_ASSERT_ERROR( eval<DenseVector<Scalar>>(p.lambda * rel_lower_bound).all( CompareOp::LE, zeroVector ), 
                       "Complementarity for lambda does not hold");

    SCAI_ASSERT_ERROR( eval<DenseVector<Scalar>>(p.lambda * rel_lower_bound).all( CompareOp::GE, zeroVector ),
                       "Complementarity for lambda does not hold");

    SCAI_ASSERT_ERROR( eval<DenseVector<Scalar>>( p.mu * rel_upper_bound ).all( CompareOp::LE, zeroVector), 
                       "Complementarity for mu does not hold.");

    SCAI_ASSERT_ERROR( eval<DenseVector<Scalar>>( p.mu * rel_upper_bound ).all( CompareOp::GE, zeroVector), 
                       "Complementarity for mu does not hold.");

    return p;
}

namespace detail
{
    template <typename Scalar>
    scai::lama::DenseVector<Scalar> randomBiasedVector(scai::IndexType size, double bias)
    {
        SCAI_ASSERT_LE(bias, 1.0, "bias must be less or equal to 1.0");
        SCAI_ASSERT_GE(bias, 0.0, "bias must be non-negative");

        std::cout << "randomBiasedVectr: " << bias << std::endl;

        using namespace scai::lama;
        const scai::dmemo::DistributionPtr dist(new scai::dmemo::NoDistribution(size));
        
        SparseVector<Scalar> v;   // 
        v.setSparseRandom( dist, 1, bias, 2 );
        v -= 1;
        return convert<DenseVector<Scalar>>( v );
    }
    
    template <typename Scalar>
    scai::lama::DenseVector<Scalar> randomBiasedNonNegativeVector(scai::IndexType size, double bias)
    {
        using namespace scai::lama;
        DenseVector<Scalar> v = randomBiasedVector<Scalar>(size, bias);
        v = abs( v );
        return v;
    }
}
