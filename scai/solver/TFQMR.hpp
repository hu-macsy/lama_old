/**
 * @file TFQMR.hpp
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
 * @brief Definion of solver class TFQMR( Transpose Free Quasi Minimal Residual)
 * @author David Schissler
 * @date 13.05.2015
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
 * @brief The class TFQMR represents a IterativeSolver which uses the krylov subspace Transpose Free
 *        Quasi Minimal Residual (TFQMR) method to solve a system of linear equations iteratively.
 *
 * Remark:
 * The scalars in the algorithm are set to zero if they are smaller then machine
 * precision (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.

 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT TFQMR:

    public IterativeSolver<ValueType>,
    public _Solver::Register<TFQMR<ValueType> >

{
public:
    /**
    * @brief Creates a TFQMR solver with a given ID.
    *
    * @param id The ID for the solver.
    */
    TFQMR( const std::string& id );
    /**
    * @brief Create a TFQMR solver with a given ID and a given logger.
    *
    * @param id        The ID of the solver.
    * @param logger    The logger which shall be used by the solver
    */
    TFQMR( const std::string& id, LoggerPtr logger );

    /**
    * @brief Copy constructor that copies the status independent solver information
    */
    TFQMR( const TFQMR& other );

    virtual ~TFQMR();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
    * @brief Copies the status independent solver informations to create a new instance of the same
    * type
    *
    * @return shared pointer of the copied solver
    */
    virtual TFQMR<ValueType>* copy();

    struct TFQMRRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        std::unique_ptr<lama::Vector<ValueType>> mVecD;
        std::unique_ptr<lama::Vector<ValueType>> mInitialR;
        std::unique_ptr<lama::Vector<ValueType>> mVecVEven;
        std::unique_ptr<lama::Vector<ValueType>> mVecVOdd;
        std::unique_ptr<lama::Vector<ValueType>> mVecVT;
        std::unique_ptr<lama::Vector<ValueType>> mVecW;
        std::unique_ptr<lama::Vector<ValueType>> mVecZ;

        RealType<ValueType> mEps;
        ValueType mAlpha;
        ValueType mBeta;
        ValueType mC;
        ValueType mEta;
        ValueType mTheta;
        ValueType mTau;
        ValueType mRhoNew;
        ValueType mRhoOld;
    };
    /**
    * @brief Returns the complete configuration of the derived class
    */
    virtual TFQMRRuntime& getRuntime();
    /**
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
    * @brief Returns the complete const configuration of the derived class
    */
    virtual const TFQMRRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    TFQMRRuntime mTFQMRRuntime;
    /**
     * @brief Performs one TFQMR iteration based on Matrix/Vector operations.
     * iterationOdd() and iterationEven() is some update for iterate() based on
     * the number of iterations (even, odd).
     */
    virtual void iterate();
    void iterationOdd();
    void iterationEven();

    using IterativeSolver<ValueType>::mPreconditioner;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
