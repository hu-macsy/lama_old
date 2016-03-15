/**
 * @file TestSolverMatrices.hpp
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
 * @brief TestSolverMatrices.hpp
 * @author David Schissler
 * @date 29.02.2016
 */

#pragma once

#include <memory>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/common/unique_ptr.hpp>
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
     * CGNR, GMRES, DefaultJacobi!!  
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
    static scai::lama::CSRSparseMatrix<ValueType> hilbertMatrix(const IndexType dim);

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
    scai::lama::MatrixCreator<ValueType>::buildPoisson2D(matrix,stencilType,dim1,dim2);
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::symmNotPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/bfwb398.mtx";

    scai::lama::CSRSparseMatrix<double> matrix(inputFile);
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::symmInDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/bcsstk22.mtx";

    scai::lama::CSRSparseMatrix<double> matrix(inputFile);
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::notSymmPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/p1femdiffusionconvection.mtx";

    scai::lama::CSRSparseMatrix<double> matrix(inputFile);
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::notSymmNotPosDefMatrix()
{
    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/west0381.mtx";

    scai::lama::CSRSparseMatrix<double> matrix(inputFile);
    return matrix;
}


template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::hilbertMatrix(const IndexType dim)
{
    ValueType *values = new ValueType[dim*dim];


    for(IndexType i=0;i<dim;++i){
        for(IndexType j=0;j<dim;++j){
            values[i*dim+j]= 1.0/(i+j+1.0);
        }
    }

    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData(dim,dim,values);
    return matrix;
}


template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSolverMatrices::complexSymmNotHermitIndefMatrix()
{
    bool isComplex;
    switch ( scai::common::TypeTraits<ValueType>::stype )
    {
        case scai::common::scalar::DOUBLE_COMPLEX:
        case scai::common::scalar::COMPLEX:
        case scai::common::scalar::LONG_DOUBLE_COMPLEX :

            isComplex = true;
            break;

        default:
            isComplex = false;
    }

        if ( isComplex == false )
    {
        COMMON_THROWEXCEPTION( "ValueType is not complex")
    }

    std::string prefix = scai::test::Configuration::getPath();
    std::string inputFile = prefix + "/p1femhelmholtz.mtx";
    scai::lama::CSRSparseMatrix<double> matrix(inputFile);
    return matrix;
}