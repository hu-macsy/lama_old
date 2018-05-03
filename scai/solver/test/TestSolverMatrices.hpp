/**
 * @file TestSolverMatrices.hpp
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
 * @brief TestSolverMatrices.hpp
 * @author David Schissler
 * @date 29.02.2016
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/common/TypeTraits.hpp>
/**
 * @brief The class TestSolverMatrices provides test matrices that
 *        are needed by the solver tests.
 *
 * The class TestSolverMatrices summaries test matrices that are needed
 * by the solver tests. TestSolverMatrices consists only of static methods
 * that return copies to Sparsematrix.
 */

class TestSolverMatrices
{
public:

    /**
     * @brief Creates a sparse, symmetric positive-definite matrix in CSR format.
     *
     * Build a sparse matrix representing the discretization of the Laplacian operator
     * on a two-dimensional structured grid. This matrix is real, sparse, symmetric and
     * positive definite.
     * It can be used to test convergence for all solvers.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] stencilType must be 5, or 9
     *  @param[in] dim1, dim2 are the sizes of the two-dimensional grid
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> buildPoisson2D(
        const IndexType stencilType,
        const IndexType dim1,
        const IndexType dim2 );

    /**
     * @brief Creates a sparse, symmetric and NOT positive-definite matrix in CSR format.
     *
     * This matrix bfwb398.mtx is real, sparse, symmetric and NOT positive definite.
     * It can be used to test convergence for the following krylow methods:
     * BiCG, CGS (slow), GMRES, MINRES, QMR
     *
     * Doesn't work for: CG, BiCGstab, TFQMR
     *
     *  @param[out] matrix in CSR format
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> symmNotPosDefMatrix();

    /**
     * @brief Creates a sparse, symmetric and indefinite matrix in CSR format.
     *
     * This matrix bcsstk22.mtx is real, sparse, symmetric, indefinite and not diagonal domainant.
     * The moduli of the entries are of order 10^6.
     *
     * Currently (14.03.2016), this leads to NaN-entries for some iterative solvers like
     * CGNR, GMRES, Jacobi!!
     *
     * It can be used to test convergence for the following krylow methods:
     * BiCG, BiCGstab, CG, CGS, MINRES, QMR
     *
     *
     *  @param[out] matrix in CSR format
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> symmInDefMatrix();

    /**
     * @brief Creates a sparse undsymmetric and positive definite matrix in CSR format.
     *
     *
     * This matrix is a unsymmetric and postive definite matrix. Since positive-
     * definiteness is defined on symmetric square matrices we say a n-dim square
     * matrix A is postive definite if x^t * A * x > 0 for all n-dim vectors x, x!=0,
     * where x^t denotes the transposition of the vector x. It is the system matrix of a
     * P1FEM discretatization of a stationary diffusion convection equation on the unit
     * square.
     * It can be used to test convergence for the following methods: splitting methods
     * Doesn't work for: krylow methods
     *
     *  @param[out] matrix in CSR format
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> notSymmPosDefMatrix();

    /**
     * @brief Creates a sparse, NOT symmetric and NOT positive-definite matrix in CSR format.
     *
     * This matrix west0381.mtx is real, sparse, unsymmetric and NOT positive definite.
     * It can be used to test convergence for the following krylow methods:
     * QMR (extremely slow)
     *
     * Doesn't work for: CG, BiCG, BiCGstab, CGS, GMRES, TFQMR
     *
     *  @param[out] matrix in CSR format
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> notSymmNotPosDefMatrix();

    /**
     * @brief Creates a dim x dim matrix in CSR format.
     *
     * This matrix is real, dense, symmetric and positive-definite. The condition
     * number grows exponentially with the dimension of the matrix. Hence,
     * krylow methods converges slowly.
     * It can be used to test preconditioners for (all solvers) krylow methods: CG, BiCG,
     * BiCGstab, CGS, GMRES, MINRES, QMR, TFQMR
     *
     * Note:
     * Currently (29.02.2016), the solvers mentioned above are limited by the
     * machine precision mEps because a small residual norm doesn't imply a small
     * error, i.e. Norm(exactSolution-approximation) >> Norm(residual).
     * Example:
     * CG solver, dim=4
     *      runtime.mEps = std::numeric_limits<double>::epsilon() * 3;
     * leads to an error while
     *      runtime.mEps = std::numeric_limits<long double>::epsilon() * 3;
     * works fine.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] dim is the size of the quadratic matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> hilbertMatrix( const IndexType dim );

    /**
     * @ Creates a complex valued matrix in CSR format.
     *
     * This matrix is complex, sparse, symmetric, not hermitian and indefinite. It is
     * the system matrix of the P1FEM discretization of some helmholtz equation on the
     * unit square equipped with some obstacles.
     * It can be used to test convergence of some krylow methods which don't assume a
     * hermitian matrix like: GMRES (complex not supported yet), BiCG, BiCGstab
     *
     *
     *  @param[out] matrix in CSR format
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> complexSymmNotHermitIndefMatrix();
};

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::buildPoisson2D(
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2 )
{
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    scai::lama::MatrixCreator<ValueType>::buildPoisson2D( matrix, stencilType, dim1, dim2 );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::symmNotPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/bfwb398.mtx";
    scai::lama::CSRSparseMatrix<double> matrix( inputFile );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::symmInDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/bcsstk22.mtx";
    scai::lama::CSRSparseMatrix<double> matrix( inputFile );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::notSymmPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/p1femdiffusionconvection.mtx";
    scai::lama::CSRSparseMatrix<double> matrix( inputFile );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::notSymmNotPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/west0381.mtx";
    scai::lama::CSRSparseMatrix<double> matrix( inputFile );
    return matrix;
}


template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::hilbertMatrix( const IndexType dim )
{
    ValueType* values = new ValueType[dim * dim];

    for ( IndexType i = 0; i < dim; ++i )
    {
        for ( IndexType j = 0; j < dim; ++j )
        {
            values[i * dim + j] = 1.0 / ( i + j + 1.0 );
        }
    }

    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( dim, dim, values );
    return matrix;
}


template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::complexSymmNotHermitIndefMatrix()
{
    bool isComplex;

    switch ( scai::common::TypeTraits<ValueType>::stype )
    {
        case scai::common::ScalarType::DOUBLE_COMPLEX:
        case scai::common::ScalarType::COMPLEX:
        case scai::common::ScalarType::LONG_DOUBLE_COMPLEX :
            isComplex = true;
            break;

        default:
            isComplex = false;
    }

    if ( isComplex == false )
    {
        COMMON_THROWEXCEPTION( "ValueType is not complex" )
    }

    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/p1femhelmholtz.mtx";
    scai::lama::CSRSparseMatrix<double> matrix( inputFile );
    return matrix;
}
