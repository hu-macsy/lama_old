/**
 * @file LAMAInputSetComplexityVisitor.hpp
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
 * @brief LAMAInputSetComplexityVisitor.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 01.10.2011
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
