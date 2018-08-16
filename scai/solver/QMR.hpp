/**
 * @file QMR.hpp
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
 * @brief QMR.hpp
 * @author lschubert
 * @date 06.08.2013
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

// scai internal libraries
#include <scai/lama/_Vector.hpp>
#include <scai/lama/Scalar.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The class QMR represents a IterativeSolver which uses the krylov subspace stabilized Quasi Minimal Residual method
 *        to solve a system of linear equations iteratively.
 *
 * This solver needs for solving a matrix A also the tranposed matrix - vector multiplication.
 * For efficiency it is often better to build the transposed matrix explicitly. Therefore
 * it is better to encapsulate a matrix A in an object AwithT.
 *
 * Remark:
 * The scalars in the algorithm are set to zero if they are smaller than machine precision
 * (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT QMR:

    public IterativeSolver<ValueType>,
    public _Solver::Register<QMR<ValueType> >

{
public:
    /**
     * @brief Creates a BiCG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    QMR( const std::string& id );

    /**
     * @brief Create a BiCG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    QMR( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    QMR( const QMR& other );

    virtual ~QMR();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual QMR<ValueType>* copy();

    struct QMRRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        lama::MatrixPtr<ValueType> mConjTransposeA;

        std::unique_ptr<lama::Vector<ValueType>> mInitialRes;
        std::unique_ptr<lama::Vector<ValueType>> mVecV;
        std::unique_ptr<lama::Vector<ValueType>> mVecW;
        std::unique_ptr<lama::Vector<ValueType>> mVecY;      /*preconditioning 1*/
        std::unique_ptr<lama::Vector<ValueType>> mVecZ;

        std::unique_ptr<lama::Vector<ValueType>> mVecWT;
        std::unique_ptr<lama::Vector<ValueType>> mVecVT;
        std::unique_ptr<lama::Vector<ValueType>> mVecYT;
        std::unique_ptr<lama::Vector<ValueType>> mVecZT;
        std::unique_ptr<lama::Vector<ValueType>> mVecP;
        std::unique_ptr<lama::Vector<ValueType>> mVecQ;
        std::unique_ptr<lama::Vector<ValueType>> mVecPT;
        std::unique_ptr<lama::Vector<ValueType>> mVecS;
        std::unique_ptr<lama::Vector<ValueType>> mVecD;

        ValueType mGamma;
        ValueType mTheta;
        ValueType mPsi;
        ValueType mRho;
        ValueType mEpsilon;
        ValueType mEta;

        ValueType mEps;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual QMRRuntime& getRuntime();
    /**
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const QMRRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    virtual void iterate();

    using IterativeSolver<ValueType>::mPreconditioner;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    QMRRuntime    mQMRRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
