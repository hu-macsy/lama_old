/**
 * @file Jacobi.hpp
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
 * @brief Jacobi.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/OmegaSolver.hpp>

// local library
#include <scai/lama/matrix/SparseMatrix.hpp>

namespace scai
{

namespace solver
{

class COMMON_DLL_IMPORTEXPORT Jacobi:
    public OmegaSolver,
    public Solver::Register<Jacobi>
{
public:
    Jacobi( const std::string& id );
    Jacobi( const std::string& id, lama::Scalar omega );
    Jacobi( const std::string& id, LoggerPtr logger );
    Jacobi( const std::string& id, lama::Scalar omega, LoggerPtr logger );
    Jacobi( const Jacobi& other );

    virtual ~Jacobi();

    virtual void initialize( const lama::_Matrix& coefficients );

    virtual void solveInit( lama::_Vector& solution, const lama::_Vector& rhs );

    virtual void solveFinalize();

    void iterate();

    struct JacobiRuntime: OmegaSolverRuntime
    {
        JacobiRuntime();
        virtual ~JacobiRuntime();

        //TODO: HArray?
        lama::_VectorPtr mOldSolution;
        SolutionProxy mProxyOldSolution;
        common::shared_ptr<hmemo::_HArray> mDiagonal;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual JacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const JacobiRuntime& getConstRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    JacobiRuntime mJacobiRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

private:
    template<typename ValueType>
    void iterateTyped( const lama::SparseMatrix<ValueType>& );

    template<typename ValueType>
    void iterateSync(
        hmemo::HArray<ValueType>& solution,
        const lama::SparseMatrix<ValueType>& coefficients,
        const hmemo::HArray<ValueType>& localOldSolution,
        hmemo::HArray<ValueType>& haloOldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega );

    template<typename ValueType>
    void iterateAsync(
        hmemo::HArray<ValueType>& solution,
        const lama::SparseMatrix<ValueType>& coefficients,
        const hmemo::HArray<ValueType>& localOldSolution,
        hmemo::HArray<ValueType>& haloOldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega );
};

} /* end namespace solver */

} /* end namespace scai */
