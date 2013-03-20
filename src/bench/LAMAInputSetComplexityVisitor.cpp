/**
 * @file LAMAInputSetComplexityVisitor.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief LAMAInputSetComplexityVisitor.cpp
 * @author jiri
 * @date 07.05.2010
 * $Id$
 */

#include <bench/LAMAInputSetComplexityVisitor.hpp>
#include <cmath>

#include <lama/matrix/SparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>

using namespace lama;

const std::string& LAMAInputSetComplexityVisitor::getGroupId( const Group group )
{
    static std::string ids[NumGroups];

    switch( group )
    {
    case CSRSAMGSpMV:
        ids[CSRSAMGSpMV] = "SparseMatrix*Vector - CSR SAMG format";
        break;
    case CSRSpMVredist:
        ids[CSRSpMVredist] = "SparseMatrix*Vector redistributed - CSR SAMG format";
        break;
    case COOSpMVredist:
        ids[COOSpMVredist] = "SparseMatrix*Vector redistributed - COO format";
        break;
    case ELLSpMVredist:
        ids[ELLSpMVredist] = "SparseMatrix*Vector redistributed - ELL format";
        break;
    case CSRYSpMV:
        ids[CSRYSpMV] = "SparseMatrix*Vector - CSR Yale format";
        break;
    case ELLSpMV:
        ids[ELLSpMV] = "SparseMatrix*Vector - ELLPACK format";
        break;
    case BCSRSpMV:
        ids[BCSRSpMV] = "SparseMatrix*Vector - BCSR format";
        break;
    case DIASpMV:
        ids[DIASpMV] = "SparseMatrix*Vector - DIA format";
        break;
    case JDSSpMV:
        ids[JDSSpMV] = "SparseMatrix*Vector - JDS format";
        break;
    case COOSpMV:
        ids[COOSpMV] = "SparseMatrix*Vector - COO format";
        break;
    case DenseMV:
        ids[DenseMV] = "DenseMatrix*Vector";
        break;
    case SpVdotDV:
        ids[SpVdotDV] = "SparseVector*Vector";
        break;
    case VdotV:
        ids[VdotV] = "Vector*Vector";
        break;
    case CSRAMGV2Jac:
        ids[CSRAMGV2Jac] = "CSR AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case ELLAMGV2Jac:
        ids[ELLAMGV2Jac] = "ELL AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case DIAAMGV2Jac:
        ids[DIAAMGV2Jac] = "DIA AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case CSRCGAMGV2Jac:
        ids[CSRCGAMGV2Jac] = "CSR CG AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case ELLCGAMGV2Jac:
        ids[ELLCGAMGV2Jac] = "ELL CG AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case DIACGAMGV2Jac:
        ids[DIACGAMGV2Jac] = "DIA CG AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case JDSCGAMGV2Jac:
        ids[JDSCGAMGV2Jac] = "JDS CG AMG V-Cycle with 2 Jacobi smoothing steps";
        break;
    case CSRSIMPLEAMGCG:
        ids[CSRSIMPLEAMGCG] = "CSR SimpleAMG preconditioned CG";
        break;
    case ELLSIMPLEAMGCG:
        ids[ELLSIMPLEAMGCG] = "ELL SimpleAMG preconditioned CG";
        break;
    case DIASIMPLEAMGCG:
        ids[DIASIMPLEAMGCG] = "DIA SimpleAMG preconditioned CG";
        break;
    case JDSSIMPLEAMGCG:
        ids[JDSSIMPLEAMGCG] = "JDS SimpleAMG preconditioned CG";
        break;
    case ELLCSRSIMPLEAMGCG:
        ids[ELLCSRSIMPLEAMGCG] = "ELL CSR SimpleAMG preconditioned CG";
        break;
    case CSRSIMPLEAMGCOMPONENTS:
        ids[CSRSIMPLEAMGCOMPONENTS] = "CSR SimpleAMG Components";
        break;
    case ELLSIMPLEAMGCOMPONENTS:
        ids[ELLSIMPLEAMGCOMPONENTS] = "ELL SimpleAMG Components";
        break;
    case DIASIMPLEAMGCOMPONENTS:
        ids[DIASIMPLEAMGCOMPONENTS] = "DIA SimpleAMG Components";
        break;
    case JDSSIMPLEAMGCOMPONENTS:
        ids[JDSSIMPLEAMGCOMPONENTS] = "JDS SimpleAMG Components";
        break;
    case CSRSIMPLEAMGSETUP:
        ids[CSRSIMPLEAMGSETUP] = "CSR SimpleAMG Setup";
        break;
    case ELLSIMPLEAMGSETUP:
        ids[ELLSIMPLEAMGSETUP] = "ELL SimpleAMG Setup";
        break;
    case DIASIMPLEAMGSETUP:
        ids[DIASIMPLEAMGSETUP] = "DIA SimpleAMG Setup";
        break;
    case JDSSIMPLEAMGSETUP:
        ids[DIASIMPLEAMGSETUP] = "JDS SimpleAMG Setup";
        break;
    case InverseSolver:
        ids[InverseSolver] = "Inverse solver";
        break;
    case LUSolver:
        ids[LUSolver] = "LU solver";
        break;
    case PLUSolver:
        ids[PLUSolver] = "PLU solver";
        break;
    case CSRJacobi:
        ids[CSRJacobi] = "CSR Jacobi";
        break;
    case ELLJacobi:
        ids[ELLJacobi] = "ELL Jacobi";
        break;
    case DIAJacobi:
        ids[DIAJacobi] = "DIA Jacobi";
        break;
    case PCSRJacobi:
        ids[PCSRJacobi] = "CSR MPI Jacobi";
        break;
    case PELLJacobi:
        ids[PELLJacobi] = "ELL MPI Jacobi";
        break;
    case PJDSJacobi:
        ids[PJDSJacobi] = "JDS MPI Jacobi";
        break;
    case MatrixBasedCSRJacobi:
        ids[MatrixBasedCSRJacobi] = "matrix based CSR Jacobi";
        break;
    case CSRPCG:
        ids[CSRPCG] = "CSR PCG";
        break;
    case ELLPCG:
        ids[ELLPCG] = "ELL PCG";
        break;
    case JDSPCG:
        ids[JDSPCG] = "JDS PCG";
        break;
    case COOPCG:
        ids[COOPCG] = "COO PCG";
        break;
    case DIAPCG:
        ids[DIAPCG] = "DIA PCG";
        break;
    case CSRPGMRES:
        ids[CSRPGMRES] = "CSR PGMRES";
        break;
    case ELLPGMRES:
        ids[ELLPGMRES] = "ELL PGMRES";
        break;
    case DIAPGMRES:
        ids[DIAPGMRES] = "DIA PGMRES";
        break;
    case SAXPY:
        ids[SAXPY] = "SAXPY";
        break;
    case VecAssignVecMinusVec:
        ids[VecAssignVecMinusVec] = "DenseVector = DenseVector - DenseVector";
        break;
    case VecAssignAdd:
        ids[VecAssignAdd] = "DenseVector += DenseVector";
        break;
    case VecAssignScale:
        ids[VecAssignScale] = "DenseVector *= scalar";
        break;
    case VecDotproduct:
        ids[VecDotproduct] = "DenseVector * DenseVector";
        break;
    case VecAssignCSRMultExpr:
    {
        ids[VecAssignCSRMultExpr] = "Expression: DenseVector = CSR Matrix * DenseVector";
    }
    break;
    case VecAssignAddCSRMultVec:
    {
        ids[VecAssignAddCSRMultVec] = "Expression: DenseVector += CSR Matrix * DenseVector";
    }
    break;
    case VecAssignVecMinusCSRMultVecExpr:
    {
        ids[VecAssignVecMinusCSRMultVecExpr] = "Expression: DenseVector = DenseVector - CSR Matrix * DenseVector";
    }
    break;
    case VecCtorVecMinusCSRMultVecExpr:
    {
        ids[VecCtorVecMinusCSRMultVecExpr] = "Expression: DenseVector[DenseVector - CSR Matrix * DenseVector]";
    }
    break;
    case GEMM:
        ids[GEMM] = "GEMM";
        break;
    case GEMV:
        ids[GEMV] = "GEMV";
        break;
    case AXPY:
        ids[AXPY] = "AXPY";
        break;
    case TRSV:
        ids[TRSV] = "TRSV";
        break;
    case TRSM:
        ids[TRSM] = "TRSM";
        break;
    case GER:
        ids[GER] = "GER";
        break;
    case SCAL:
        ids[SCAL] = "SCAL";
        break;
    case IXAMAX:
        ids[IXAMAX] = "IXAMAX";
        break;
    case ConvertCSR2JDS:
        ids[ConvertCSR2JDS] = "CSR 2 JDS Conversion";
        break;
    case CSRPMetaSolver:
        ids[CSRPMetaSolver] = "PMetaSolver - CSR format";
        break;
    case ELLPMetaSolver:
        ids[ELLPMetaSolver] = "PMetaSolver - ELL format";
        break;
    case DIAPMetaSolver:
        ids[DIAPMetaSolver] = "PMetaSolver - DIA format";
        break;
    default:
        //This will never happen, it is only here because to have a
        //sense full default branch.
        ids[NumGroups - 1] = "Unknown Group";
        break;
    }

    return ids[group];
}

namespace
{

template<typename IndexType>
void calcSpMVComplexity(
    const IndexType numRows,
    const IndexType numValues,
    CounterType& numFlops,
    CounterType& processedBytesFloat,
    CounterType& processedBytesDouble )
{
    //every matrix element is multiplied once
    //and all products of one row of the matrix are added together
    //so we have nna float multiplies and nna - nnu float adds.
    numFlops = 2 * numValues - numRows;
    //The whole matrix need to be accessed once
    //for each row of the matrix two values of the index array need to be loaded
    CounterType processedBytes = ( numValues + 2 * numRows ) * sizeof(IndexType);
    //we need to load the values of the matrix once
    CounterType numProcessedValues = numValues;
    //we need to load one value of the input vector for each element of
    //the matrix (we ignore the cache)
    numProcessedValues += numRows;
    //we need to write each value of the output vector once
    numProcessedValues += numRows;
    processedBytesFloat = processedBytes;
    processedBytesFloat += numProcessedValues * sizeof(float);
    processedBytesDouble = processedBytes;
    processedBytesDouble += numProcessedValues * sizeof(double);
}

}

void LAMAInputSetComplexityVisitor::getMVComplexity(
    const Matrix& matrix,
    CounterType& numFlops,
    CounterType& processedBytesFloat,
    CounterType& processedBytesDouble )
{
    IndexType numRows = matrix.getNumRows(); // used for all formats

    if( matrix.getMatrixKind() == Matrix::DENSE )
    {
        IndexType numCols = matrix.getNumColumns();
        numFlops = numRows * ( 2 * numCols - 1 );
        CounterType processedValues = 2 * numRows * numCols + numRows;
        processedBytesFloat = processedValues * sizeof(float);
        processedBytesDouble = processedValues * sizeof(double);
        return;
    }

    const _SparseMatrix* sparseM = dynamic_cast<const _SparseMatrix*>( &matrix );

    if( sparseM )
    {
        const _MatrixStorage& localData = sparseM->getLocalStorage();
        const _MatrixStorage& haloData = sparseM->getHaloStorage();
        IndexType numLocalValues = localData.getNumValues() + haloData.getNumValues();
        const Communicator& comm = matrix.getDistribution().getCommunicator();
        IndexType numValues = comm.sum( numLocalValues );
        //every matrix element is multiplied once
        //and all products of one row of the matrix are added together
        //so we have numValues float multiplies and numValues - numRows float adds.
        numFlops = 2 * numValues - numRows;
        //The whole matrix need to be accessed once
        //for each row of the matrix two values of the index array need to be loaded
        CounterType processedBytes = ( numValues + numRows ) * sizeof(IndexType);
        //we need to load the values of the matrix once
        CounterType numProcessedValues = numValues;
        //we need to load one value of the input vector for each element of
        //the matrix (we ignore the cache)
        numProcessedValues += numValues;
        //we need to write each value of the output vector once
        numProcessedValues += numRows;
        processedBytesFloat = processedBytes;
        processedBytesFloat += numProcessedValues * sizeof(float);
        processedBytesDouble = processedBytes;
        processedBytesDouble += numProcessedValues * sizeof(double);
        return;
    }

    LAMA_THROWEXCEPTION( "unrecognized matrix " << matrix << ": no MV complexity available" );
}

void LAMAInputSetComplexityVisitor::getCGComplexity(
    const Matrix& matrix,
    CounterType& numFlops,
    CounterType& processedBytesFloat,
    CounterType& processedBytesDouble )
{
    // init to 0
    numFlops = 0;
    processedBytesFloat = 0;
    processedBytesDouble = 0;

    // Total work is:
    // 1 - residual computation
    // 2 - vec = vec
    // 3 - (vec,vec)  (2 times)
    // 4 - vec = mat * Vec
    // 5 - vec = vec +- alpha*vec (3 times)
    CounterType bytesFloat = 0;
    CounterType bytesDouble = 0;
    const IndexType numRows = matrix.getNumRows();
    //Complexity of residual calculation SpMV (1+4)
    getMVComplexity( matrix, numFlops, bytesFloat, bytesDouble );
    //Calculation of the residual is v = y - A*x not SpMV (v = A*x)
    numFlops += ( 2 * numFlops + numRows );
    processedBytesFloat += ( 2 * bytesFloat + numRows * sizeof(float) );
    processedBytesDouble += ( 2 * bytesDouble + numRows * sizeof(double) );
    //Complexity of Vector Assignment vec=vec (2)
    processedBytesFloat += ( 2 * numRows * sizeof(float) );
    processedBytesDouble += ( 2 * numRows * sizeof(double) );
    //Complexity of 3 Vector Assignment SAXPY(5)
    numFlops += 3 * ( 2 * numRows );
    processedBytesFloat += 3 * ( ( 3 * numRows + 1 ) * sizeof(float) );
    processedBytesDouble += 3 * ( ( 3 * numRows + 1 ) * sizeof(double) );
    //Complexity of 2 (,) product Assignment (3)
    numFlops += 2 * ( 2 * numRows );
    processedBytesFloat += 2 * ( ( 2 * numRows + 1 ) * sizeof(float) );
    processedBytesDouble += 2 * ( ( 2 * numRows + 1 ) * sizeof(double) );
}

void LAMAInputSetComplexityVisitor::getGMRESComplexity(
    const Matrix&,
    CounterType& numFlops,
    CounterType& processedBytesFloat,
    CounterType& processedBytesDouble )
{
    // init to 0
    numFlops = 0;
    processedBytesFloat = 0;
    processedBytesDouble = 0;

    //TODO needs to be implemented as soon as algorithm remains unchanged!!!
}

void LAMAInputSetComplexityVisitor::getDefaultJacobiComplexity(
    unsigned int numIterations,
    const Matrix& matrix,
    CounterType& numFlops,
    CounterType& processedBytesFloat,
    CounterType& processedBytesDouble )
{
    // init to 0
    numFlops = 0;
    processedBytesFloat = 0;
    processedBytesDouble = 0;

    // Total work is:
    // 1 - compute d^-1 * rhs (only done once)
    // 2 - vec = vec- mat * vec (done numIterations times)
    // 3 - vec = w*vec + (1+w)*vec
    CounterType flops = 0;
    CounterType bytesFloat = 0;
    CounterType bytesDouble = 0;
    const IndexType numRows = matrix.getNumRows();

    //Complexity of residual calculation SpMV (1+2)
    getMVComplexity( matrix, flops, bytesFloat, bytesDouble );

    // 1 - compute d^-1 * rhs
    numFlops += numRows;
    processedBytesFloat += ( 2 * numRows * sizeof(float) );
    processedBytesDouble += ( 2 * numRows * sizeof(double) );

    // 2 - vec = vec- mat * vec (done numIterations times)
    numFlops += numIterations * ( numRows + flops );
    processedBytesFloat += numIterations * ( bytesFloat + numRows * sizeof(float) );
    processedBytesDouble += numIterations * ( bytesDouble + numRows * sizeof(double) );

    // 3 - vec = w*vec + (1+w)*vec
    numFlops += numIterations * ( 3 * numRows );
    processedBytesFloat += numIterations * ( 2 * numRows * sizeof(float) );
    processedBytesFloat += numIterations * ( 2 * numRows * sizeof(double) );
}
