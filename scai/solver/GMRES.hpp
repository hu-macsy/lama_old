/**
 * @file GMRES.hpp
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
class COMMON_DLL_IMPORTEXPORT GMRES:
    public IterativeSolver,
    public Solver::Register<GMRES>
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

    virtual void initialize( const lama::Matrix& coefficients );

    void setKrylovDim( unsigned int krylovDim );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct GMRESRuntime: IterativeSolverRuntime
    {
        GMRESRuntime();
        virtual ~GMRESRuntime();

        // arrays to store rotations
        std::unique_ptr<double[]> mCC;
        std::unique_ptr<double[]> mSS;

        // array for Hessenberg equation
        // H*y=g
        std::unique_ptr<double[]> mG;
        std::unique_ptr<double[]> mY;

        // Hessenberg matrix
        // mH:  Upper triangular (columnwise)
        // mHd: diagonal band h(i+1,i)
        std::unique_ptr<double[]> mH;
        std::unique_ptr<double[]> mHd;

        // krylov space
        std::vector<lama::Vector*>* mV;

        // temp-arrays
        lama::Vector* mW;
        lama::Vector* mT;

        // remember starting solution
        // only needed if x is modified within krylov loop
        lama::Vector* mX0;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual GMRESRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const GMRESRuntime& getConstRuntime() const;

    double getAverageIterationTime() const;
    double getAveragePreconditionerTime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    virtual void iterate();

    GMRESRuntime mGMRESRuntime;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    void updateX( unsigned int i );

    // krylov dimension
    unsigned int mKrylovDim;

    double totalIterationTime;
    double totalPreconditionerTime;
};

} /* end namespace solver */

} /* end namespace scai */
