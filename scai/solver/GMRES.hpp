/**
 * @file GMRES.hpp
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
 * @brief GMRES.hpp
 * @author Malte Förster
 * @date 10.04.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

#include <memory>

namespace scai
{

namespace solver
{

/**
 * @brief The class GMRES represents a IterativeSolver which uses the krylov subspace GMRES method
 *        to solve a system of linear equations iteratively.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT GMRES:

    public IterativeSolver<ValueType>,
    public _Solver::Register<GMRES<ValueType> >
{
public:

    /**
     * @brief Creates a GMRES solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    GMRES( const std::string& id );

    /**
     * @brief Create a GMRES solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    GMRES( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    GMRES( const GMRES& other );

    virtual ~GMRES();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    void setKrylovDim( IndexType krylovDim );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual GMRES<ValueType>* copy();

    /** Runtime data for the GMRES solver.
     *
     *  Most data strucures are arrays of size krylovDim + 1, all are implemented
     *  via unique_ptr to guarantee correct delete in the destructor.
     */
    struct GMRESRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        // arrays to store rotations
        std::vector<ValueType> mCC;
        std::vector<ValueType> mSS;

        // array for Hessenberg equation
        // H*y=g
        std::vector<ValueType> mG;
        std::vector<ValueType> mY;

        // Hessenberg matrix
        // mH:  Upper triangular (columnwise)
        // mHd: diagonal band h(i+1,i)
        std::vector<ValueType> mH;
        std::vector<ValueType> mHd;

        // krylov space, vector of (unique) vector pointers
        std::vector<std::unique_ptr<lama::Vector<ValueType>>> mV;

        // temp-arrays
        std::unique_ptr<lama::Vector<ValueType>> mW;
        std::unique_ptr<lama::Vector<ValueType>> mT;

        // remember starting solution
        // only needed if x is modified within krylov loop
        std::unique_ptr<lama::Vector<ValueType>> mX0;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual GMRESRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const GMRESRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    virtual void iterate();

    GMRESRuntime mGMRESRuntime;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    using IterativeSolver<ValueType>::mPreconditioner;

private:

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    void updateX( IndexType i );

    // krylov dimension
    IndexType mKrylovDim;
};

} /* end namespace solver */

} /* end namespace scai */
