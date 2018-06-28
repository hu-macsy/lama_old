/**
 * @file ipbcls/NewtonStepCG.hpp
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
 * @brief ToDo: Missing description in ./ipbcls/NewtonStepCG.hpp
 * @author Thomas.Brandes@scai.fraunhofer.de 2018-03-02
 * @date 16.03.2015
 */
#include <scai/lama/Vector.hpp>
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/solver/Solver.hpp>

/**
 * NewtonStepCG computes search directions for the truncated Newton's method,
 * as used in optimization.
 *
 * TODO: Docs!
 */
template<typename ValueType>
class NewtonStepCG
{
public:

    NewtonStepCG( const scai::lama::Matrix<ValueType>* hessian);

    ValueType getForcingTerm() const;
    void setForcingTerm(const ValueType forcingTerm);

    scai::IndexType getMaxIterations() const;
    void setMaxIterations(scai::IndexType maxIter);

    scai::solver::SolverPtr<ValueType> getPreconditioner() const;
    void setPreconditioner(scai::solver::SolverPtr<ValueType> hessianPreconditioner);

    /**
     * Computes an appropriate search direction for the current gradient and Hessian.
     */
    scai::IndexType computeStep(scai::lama::Vector<ValueType> & dx, const scai::lama::Vector<ValueType> & gradient) const;

private:

    scai::lama::Vector<ValueType> & precondition(const scai::lama::Vector<ValueType> & r) const;

    const scai::lama::Matrix<ValueType> * mHessian;
    ValueType mForcingTerm;
    scai::IndexType mMaxIterations;
    scai::solver::SolverPtr<ValueType> mHessianPreconditioner;

    // Workspace variables, must be pointers as Vector<ValueType> is abstract

    std::unique_ptr<scai::lama::Vector<ValueType>> rPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> zPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> pPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> HpPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> rMinusGradientPtr;
};
