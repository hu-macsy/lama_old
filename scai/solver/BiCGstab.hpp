/**
 * @file BiCGstab.hpp
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
 * @brief BiCGstab.hpp
 * @author Lauretta Schubert
 * @date 06.08.2013
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>
#include <scai/lama/Vector.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The class BiCGstab represents a IterativeSolver which uses the krylov subspace stabilized BiCG method
 *        to solve a system of linear equations iteratively.
 *
 * Remark:
 * The scalars in the algorithm are set to zero if they are smaller than machine precision
 * (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BiCGstab:

    public IterativeSolver<ValueType>,
    public _Solver::Register<BiCGstab<ValueType> >
{
public:
    /**
     * @brief Creates a BiCG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    BiCGstab( const std::string& id );

    /**
     * @brief Create a BiCG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    BiCGstab( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    BiCGstab( const BiCGstab& other );

    virtual ~BiCGstab();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual BiCGstab<ValueType>* copy();

    /** 
     *  @brief Runtime variables for BiCGstab solver.
     *
     *  Note: due to polymorphism we use unique pointers for temporary vectors
     *        (will be DenseVector in most cases).
     */
    struct BiCGstabRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        std::unique_ptr<lama::Vector<ValueType>> mRes0;
        std::unique_ptr<lama::Vector<ValueType>> mVecV;
        std::unique_ptr<lama::Vector<ValueType>> mVecP;
        std::unique_ptr<lama::Vector<ValueType>> mVecS;
        std::unique_ptr<lama::Vector<ValueType>> mVecT;
        std::unique_ptr<lama::Vector<ValueType>> mVecPT;
        std::unique_ptr<lama::Vector<ValueType>> mVecST;
        std::unique_ptr<lama::Vector<ValueType>> mVecTT;

        RealType<ValueType> mEps;      // used in comparison
        RealType<ValueType> mResNorm;  // norm is always real
        ValueType mOmega;
        ValueType mAlpha;
        ValueType mBeta;
        ValueType mRhoOld;
        ValueType mRhoNew;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGstabRuntime& getRuntime();

    /**
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const BiCGstabRuntime& getRuntime() const;

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

    BiCGstabRuntime    mBiCGstabRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
