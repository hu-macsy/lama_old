/**
 * @file solver/examples/myJacobi/MyJacobiModule.hpp
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
 * @brief MyJacobi.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/OmegaSolver.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/DenseVector.hpp>

// logging
#include <scai/logging/Logger.hpp>

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MyJacobi:

    public scai::solver::OmegaSolver<ValueType>,
    public scai::solver::_Solver::Register<MyJacobi<ValueType> >
{
public:

    MyJacobi( const std::string& id );

    MyJacobi( const std::string& id, scai::solver::LoggerPtr logger );

    MyJacobi( const std::string& id, const ValueType omega ); 

    MyJacobi( const std::string& id, const ValueType omega, scai::solver::LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    MyJacobi( const MyJacobi<ValueType>& other );

    virtual ~MyJacobi();

    /**
     * @brief Initializes the solver by calculating D^(-1) * C from A = D + C.
     *
     * This specialized initialization for the Jacobi solver does the
     * one-time calculation of D^(-1)*C where D is the diagonal Matrix
     * of A and C is the rest (A = D + C).
     *
     * @param coefficients The matrix A from A*u=f
     */
    virtual void initialize( const scai::lama::Matrix<ValueType>& coefficients );

    /** 
     * @brief Override default implementation 
     */
    virtual void solveInit( scai::lama::Vector<ValueType>& solution, const scai::lama::Vector<ValueType>& rhs );

    /**
     * @brief Implementation of pure copy method with covariant return type.
     */
    virtual MyJacobi* copy();

    using typename scai::solver::IterativeSolver<ValueType>::IterativeSolverRuntime;

    struct MyJacobiRuntime:scai::solver::IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        MyJacobiRuntime();

        virtual ~MyJacobiRuntime();

        scai::lama::MatrixPtr<ValueType>   mDiagonalTimesLU;   // modified copy of A required
        scai::lama::DenseVector<ValueType> mDiagonalInverted;  // keeps 1 / D
        scai::lama::DenseVector<ValueType> mDiagonalTimesRhs;  // rhs * D, only computed once
        scai::lama::DenseVector<ValueType> mOldSolution;       // temporary, used to swap with solution
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual MyJacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const MyJacobiRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factory

    static scai::solver::SolverCreateKeyType createValue();

    // static method for create by factory

    static scai::solver::_Solver* create();

protected:

    MyJacobiRuntime mMyJacobiRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Implementation of pure methode Iterativesolver<ValueType>::iterate
     *
     * In addition to the Jacobi iteration the term D^(-1)*f from A=u*f and
     * A = D + C is evaluated before the first iteration.
     */
    virtual void iterate();
};
