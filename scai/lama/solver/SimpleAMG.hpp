/**
 * @file SimpleAMG.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief SimpleAMG.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/Solver.hpp>
#include <scai/lama/solver/IterativeSolver.hpp>

// local library
#include <scai/lama/solver/AMGSetup.hpp>

// std
#include <vector>

namespace scai
{

namespace lama
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

    virtual void initialize( const Matrix& coefficients );

    virtual void iterate();

    void setMaxLevels( unsigned int levels );

    void setMinVarsCoarseLevel( unsigned int vars );

    unsigned int getNumLevels();

    const Matrix& getGalerkin( unsigned int level );
    const Matrix& getRestriction( unsigned int level );
    const Matrix& getInterpolation( unsigned int level );

    Vector& getSolutionVector( unsigned int level );
    Vector& getRhsVector( unsigned int level );

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

} /* end namespace lama */

} /* end namespace scai */
