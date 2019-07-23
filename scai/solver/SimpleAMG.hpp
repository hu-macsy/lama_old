/**
 * @file SimpleAMG.hpp
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

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SimpleAMG:

    public IterativeSolver<ValueType>,
    public _Solver::Register<SimpleAMG<ValueType> >
{
public:

    SimpleAMG( const std::string& id );

    SimpleAMG( const std::string& id, LoggerPtr logger );

    SimpleAMG( const SimpleAMG& other );

    virtual ~SimpleAMG();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    virtual void iterate();

    void setMaxLevels( IndexType levels );

    void setMinVarsCoarseLevel( IndexType vars );

    IndexType getNumLevels();

    const lama::Matrix<ValueType>& getGalerkin( IndexType level );
    const lama::Matrix<ValueType>& getRestriction( IndexType level );
    const lama::Matrix<ValueType>& getInterpolation( IndexType level );

    lama::Vector<ValueType>& getSolutionVector( IndexType level );
    lama::Vector<ValueType>& getRhsVector( IndexType level );

    Solver<ValueType>& getSmoother( IndexType level );
    Solver<ValueType>& getCoarseLevelSolver();

    void setSmootherContext( hmemo::ContextPtr smootherContext );

    void setHostOnlyLevel( IndexType hostOnlyLevel );

    void setHostOnlyVars( IndexType hostOnlyVars );

    void setReplicatedLevel( IndexType replicatedLevel );

    void setCoarseLevelSolver( SolverPtr<ValueType> solver );

    /**
     * @brief Sets the smoother for all level
     */
    void setSmoother( SolverPtr<ValueType> solver );

    struct SimpleAMGRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        SimpleAMGRuntime();
        ~SimpleAMGRuntime();

        IndexType mCurrentLevel;  // used for recursive calls in cycle
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual SimpleAMGRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SimpleAMGRuntime& getRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same type
     *
     * @return shared pointer of the copied solver
     */
    virtual SimpleAMG<ValueType>* copy();

    double getAverageSmootherTime() const;
    double getAverageTransferTime() const;
    double getAverageResidualTime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

    /**
     *  @brief Getter method to get a const reference to the AMGSetup for queries.
     */
    const AMGSetup<ValueType>& getSetup() const;

protected:

    SimpleAMGRuntime mSimpleAMGRuntime;

    using Solver<ValueType>::mLogger;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

private:

    std::unique_ptr<AMGSetup<ValueType> > mSetup;

    void logSetupSettings();
    void logSetupInfo();
    void logSolverInfo();
    void logSetupDetails();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    void cycle();

    static void loadSetupLibs();

    /** 
     *   Help function that creates a new setup corresponding to the environment variables
     */
    static AMGSetup<ValueType>* newSetup();

    double totalSmootherTime;
    double totalTransferTime;
    double totalResidualTime;
    int totalIterations;
};

} /* end namespace solver */

} /* end namespace scai */
