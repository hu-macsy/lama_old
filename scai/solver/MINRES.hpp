/**
 * @file MINRES.hpp
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
 * @brief MINRES.hpp
 * @author David Schissler
 * @date 29.05.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{
/**
 * @brief The class MINRES represents a IterativeSolver which uses the krylov subspace Minimum Residual (MINRES)
 * method to solve a system of linear equations iteratively.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MINRES:
    public IterativeSolver<ValueType>,
    public _Solver::Register<MINRES<ValueType> >
{
public:
    /**
    * @brief Creates a MINRES solver with a given ID.
    *
    * @param id The ID for the solver.
    */
    MINRES( const std::string& id );
    /**
    * @brief Create a MINRES solver with a given ID and a given logger.
    *
    * @param id        The ID of the solver.
    * @param logger    The logger which shall be used by the solver
    */
    MINRES( const std::string& id, LoggerPtr logger );

    /**
    * @brief Copy constructor that copies the status independent solver information
    */
    MINRES( const MINRES& other );

    virtual ~MINRES();

    /** @brief Initialize the solver specific runtime data when a matrix is set.  */

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
    * @brief Implementation of pure method to create a copy of this solver 
    *        (without copying the runtime data).
    *
    * @return raw pointer of the copied solver, caller takes ownership
    */
    virtual MINRES<ValueType>* copy();

    struct MINRESRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        std::unique_ptr<lama::Vector<ValueType>> mVecV;
        std::unique_ptr<lama::Vector<ValueType>> mVecVOld;
        std::unique_ptr<lama::Vector<ValueType>> mVecVNew;
        std::unique_ptr<lama::Vector<ValueType>> mVecP;
        std::unique_ptr<lama::Vector<ValueType>> mVecPOld;
        std::unique_ptr<lama::Vector<ValueType>> mVecPNew;

        ValueType mAlpha;
        RealType<ValueType> mBetaNew;
        RealType<ValueType> mBeta;
        ValueType mC;
        ValueType mCOld;
        ValueType mCNew;
        ValueType mS;
        ValueType mSOld;
        ValueType mSNew;
        ValueType mZeta;

        RealType<ValueType> mEps;
    };
    /**
    * @brief Returns the complete configuration of the derived class
    */
    virtual MINRESRuntime& getRuntime();
    /**
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
    * @brief Returns the complete const configuration of the derived class
    */
    virtual const MINRESRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    MINRESRuntime mMINRESRuntime;
    /**
     * @brief Implementation of pure method IterativeSolver<ValueType>::iterate
     * 
     * Performs one MINRES iteration based on Matrix/Vector operations
     */
    virtual void iterate();

    void Lanczos();

    void applyGivensRotation();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;
};

} /* end namespace solver */

} /* end namespace scai */
