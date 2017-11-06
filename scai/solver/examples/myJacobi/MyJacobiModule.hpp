/**
 * @file solver/examples/myJacobi/MyJacobiModule.hpp
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

// logging
#include <scai/logging/Logger.hpp>

class COMMON_DLL_IMPORTEXPORT MyJacobi:
    public scai::solver::OmegaSolver,
    public scai::solver::Solver::Register<MyJacobi>
{
public:
    MyJacobi( const std::string& id );

    MyJacobi( const std::string& id, scai::solver::LoggerPtr logger );

    MyJacobi( const std::string& id, const scai::lama::Scalar omega ); //2nd param _Matrix.Scalar

    MyJacobi( const std::string& id, const scai::lama::Scalar omega, scai::solver::LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    MyJacobi( const MyJacobi& other );

    virtual ~MyJacobi();

    /**
     * @brief Initializes the solver by calculating D^(-1)*C from A = D + C.
     *
     * This specialized initialization for the Jacobi solver does the
     * one-time calculation of D^(-1)*C where D is the diagonal _Matrix
     * of A and C is the rest (A = D + C).
     *
     * @param coefficients The matrix A from A*u=f
     */
    virtual void initialize( const scai::lama::_Matrix& coefficients );

    virtual void solve( scai::lama::_Vector& solution, const scai::lama::_Vector& rhs );

    virtual void solveInit( scai::lama::_Vector& solution, const scai::lama::_Vector& rhs );

    virtual void solveFinalize();

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual scai::solver::SolverPtr copy();

    struct MyJacobiRuntime: OmegaSolverRuntime
    {
        MyJacobiRuntime();
        virtual ~MyJacobiRuntime();

        scai::lama::_MatrixPtr mDiagonalTimesLU;
        scai::lama::_MatrixPtr mDiagonalInverted;
        scai::lama::_VectorPtr mDiagonalTimesRhs;
        scai::lama::_VectorPtr mOldSolution;

        scai::solver::SolutionProxy mProxyOldSolution;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual MyJacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const MyJacobiRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:
    MyJacobiRuntime mMyJacobiRuntime;

    /**
     * @brief Performs one Jacobi iteration based on _Matrix/Vector operations
     *
     * In addition to the Jacobi iteration the term D^(-1)*f from A=u*f and
     * A = D + C is evaluated before the first iteration.
     */
    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    template<typename ValueType>
    void initialize( const scai::lama::_Matrix& coefficients );

    template<typename ValueType>
    void iterate();
};
