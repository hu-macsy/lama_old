/**
 * @file LAMAInputSetComplexityVisitor.hpp
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
 * @brief LAMAInputSetComplexityVisitor.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 01.10.2011
 * $Id$
 */

#pragma once

#include <scai/lama/benchmark/LAMAInputSet.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/lama/matrix/Matrix.hpp>

class LAMAInputSetComplexityVisitor
{
public:
    enum Group
    {
        CSRSAMGSpMV,
        CSRSpMVredist,
        COOSpMVredist,
        ELLSpMVredist,
        CSRYSpMV,
        ELLSpMV,
        BCSRSpMV,
        DIASpMV,
        JDSSpMV,
        COOSpMV,
        DenseMV,
        SpVdotDV,
        VdotV,
        CSRAMGV2Jac,
        ELLAMGV2Jac,
        DIAAMGV2Jac,
        CSRCGAMGV2Jac,
        ELLCGAMGV2Jac,
        DIACGAMGV2Jac,
        JDSCGAMGV2Jac,
        CSRSIMPLEAMGCG,
        ELLSIMPLEAMGCG,
        DIASIMPLEAMGCG,
        JDSSIMPLEAMGCG,
        ELLCSRSIMPLEAMGCG,
        CSRSIMPLEAMGCOMPONENTS,
        ELLSIMPLEAMGCOMPONENTS,
        DIASIMPLEAMGCOMPONENTS,
        JDSSIMPLEAMGCOMPONENTS,
        CSRSIMPLEAMGSETUP,
        ELLSIMPLEAMGSETUP,
        DIASIMPLEAMGSETUP,
        JDSSIMPLEAMGSETUP,
        InverseSolver,
        LUSolver,
        PLUSolver,
        CSRJacobi,
        ELLJacobi,
        DIAJacobi,
        PCSRJacobi,
        PELLJacobi,
        PJDSJacobi,
        MatrixBasedCSRJacobi,
        CSRPCG,
        ELLPCG,
        JDSPCG,
        DIAPCG,
        COOPCG,
        CSRPMetaSolver,
        ELLPMetaSolver,
        DIAPMetaSolver,
        CSRPGMRES,
        ELLPGMRES,
        DIAPGMRES,
        SAXPY,
        GEMM,
        GEMV,
        AXPY,
        SCAL,
        IXAMAX,
        TRSV,
        TRSM,
        GER,
        VecAssignAdd,
        VecAssignScale,
        VecDotproduct,
        VecAssignVecMinusVec,
        VecAssignCSRMultExpr,
        VecAssignAddCSRMultVec,
        VecAssignVecMinusCSRMultVecExpr,
        VecCtorVecMinusCSRMultVecExpr,
        ConvertCSR2JDS,
        NumGroups
    };

    static const std::string& getGroupId( const Group group );

    static void getMVComplexity(
        const scai::lama::Matrix& matrix,
        CounterType& numFlops,
        CounterType& processedBytesFloat,
        CounterType& processedBytesDouble );

    static void getCGComplexity(
        const scai::lama::Matrix& matrix,
        CounterType& numFlops,
        CounterType& processedBytesFloat,
        CounterType& processedBytesDouble );

    static void getGMRESComplexity(
        const scai::lama::Matrix& matrix,
        CounterType& numFlops,
        CounterType& processedBytesFloat,
        CounterType& processedBytesDouble );

    static void getDefaultJacobiComplexity(
        unsigned int numIterations,
        const scai::lama::Matrix& matrix,
        CounterType& numFlops,
        CounterType& processedBytesFloat,
        CounterType& processedBytesDouble );

private:

    LAMAInputSetComplexityVisitor();
    LAMAInputSetComplexityVisitor( const LAMAInputSetComplexityVisitor& other );
    virtual ~LAMAInputSetComplexityVisitor();

};
