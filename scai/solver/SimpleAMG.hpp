/**
 * @file SimpleAMG.hpp
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
 * @brief SimpleAMG.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

// local library
#include <scai/solver/AMGSetup.hpp>

// std
#include <vector>

namespace scai
{

namespace solver
{

class COMMON_DLL_IMPORTEXPORT SimpleAMG:
    public IterativeSolver,
    public Solver::Register<SimpleAMG>
{
public:

    SimpleAMG( const std::string& id );

    SimpleAMG( const std::string& id, LoggerPtr logger );

    SimpleAMG( const SimpleAMG& other );

    virtual ~SimpleAMG();

    virtual void initialize( const lama::Matrix& coefficients );

    virtual void iterate();

    void setMaxLevels( unsigned int levels );

    void setMinVarsCoarseLevel( unsigned int vars );

    unsigned int getNumLevels();

    const lama::Matrix& getGalerkin( unsigned int level );
    const lama::Matrix& getRestriction( unsigned int level );
    const lama::Matrix& getInterpolation( unsigned int level );

    lama::Vector& getSolutionVector( unsigned int level );
    lama::Vector& getRhsVector( unsigned int level );

    Solver& getSmoother( unsigned int level );
    Solver& getCoarseLevelSolver();

    void setSmootherContext( hmemo::ContextPtr smootherContext );

    void setHostOnlyLevel( IndexType hostOnlyLevel );

    void setHostOnlyVars( IndexType hostOnlyVars );

    void setReplicatedLevel( IndexType replicatedLevel );

    void setCoarseLevelSolver( SolverPtr solver );

    /**
     * @brief Sets the smoother for all level
     */
    void setSmoother( SolverPtr solver );

    struct SimpleAMGRuntime: IterativeSolverRuntime
    {
        SimpleAMGRuntime();
        virtual ~SimpleAMGRuntime();

        common::shared_ptr<AMGSetup> mSetup;
        unsigned int mCurrentLevel;
        void* mLibHandle;
        IndexType mHostOnlyLevel;
        IndexType mHostOnlyVars;
        IndexType mReplicatedLevel;

        SCAI_LOG_DECL_STATIC_LOGGER( logger )
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual SimpleAMGRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SimpleAMGRuntime& getConstRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    double getAverageSmootherTime() const;
    double getAverageTransferTime() const;
    double getAverageResidualTime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    SimpleAMGRuntime mSimpleAMGRuntime;

    unsigned int mMaxLevels;
    unsigned int mMinVarsCoarseLevel;
    SolverPtr mCoarseLevelSolver;
    SolverPtr mSmoother;
    hmemo::ContextPtr mSmootherContext;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

private:

    void logSetupSettings();
    void logSetupInfo();
    void logSolverInfo();
    void logSetupDetails();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    void cycle();

    void loadSetupLibs();

    double totalSmootherTime;
    double totalTransferTime;
    double totalResidualTime;
    int totalIterations;
};

} /* end namespace solver */

} /* end namespace scai */
