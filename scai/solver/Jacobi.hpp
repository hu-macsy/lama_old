/**
 * @file Jacobi.hpp
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
#include <scai/lama/Vector.hpp>
#include <scai/lama/DenseVector.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Jacobi:
    public OmegaSolver<ValueType>,
    public _Solver::Register<Jacobi<ValueType> >
{
public:

    Jacobi( const std::string& id );

    Jacobi( const std::string& id, ValueType omega );

    Jacobi( const std::string& id, LoggerPtr logger );

    Jacobi( const std::string& id, ValueType omega, LoggerPtr logger );

    Jacobi( const Jacobi& other );

    virtual ~Jacobi();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /** Implementation of pure method IterativeSolver<ValueType>::iterate */

    void iterate();

    /** Note: Jacobi can only deal with dense vectors */

    struct JacobiRuntime: OmegaSolver<ValueType>::IterativeSolverRuntime
    {
        lama::DenseVector<ValueType> mDiagonal;     //!< stores the diagonal of matrix to solve
        lama::DenseVector<ValueType> mOldSolution;  //!< temporary, reused in each iteration
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual JacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const JacobiRuntime& getRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual Jacobi<ValueType>* copy();

    static SolverCreateKeyType createValue();

    static _Solver* create();

protected:

    JacobiRuntime mJacobiRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
