/**
 * @file SimpleAMG.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of methods for the default AMG solver SimpleAMG.
 * @author Jiri Kraus
 * @date 27.10.2011
 */

// hpp
#include <scai/solver/SimpleAMG.hpp>

// local library
#include <scai/solver/SingleGridSetup.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

// internal scai libraries
#include <scai/common/Settings.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/Walltime.hpp>

#include <scai/tracing.hpp>

// std
#include <cstdlib>
#include <iomanip>

using namespace scai::hmemo;

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( SimpleAMG::logger, "Solver.IterativeSolver.SimpleAMG" )
SCAI_LOG_DEF_LOGGER( SimpleAMG::SimpleAMGRuntime::logger, "Solver.IterativeSolver.SimpleAMG.SimpleAMGRuntime" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

SimpleAMG::SimpleAMG( const std::string& id )
    : IterativeSolver( id ), mMaxLevels( 25 ), mMinVarsCoarseLevel( 100 ), mSmootherContext(
          Context::getHostPtr() )
{
    SCAI_LOG_INFO( logger, "SimpleAMG, id = " << id << " created, no logger" )
}

SimpleAMG::SimpleAMG( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger ), mMaxLevels( 25 ), mMinVarsCoarseLevel( 100 ), mSmootherContext(
          Context::getHostPtr() )
{
    SCAI_LOG_INFO( SimpleAMG::logger, "SimpleAMG, id = " << id << " created, with logger" )
}

SimpleAMG::SimpleAMG( const SimpleAMG& other )
    : IterativeSolver( other ), mMaxLevels( other.mMaxLevels ), mMinVarsCoarseLevel(
          other.mMinVarsCoarseLevel ), mSmootherContext( other.mSmootherContext )
{
}

SimpleAMG::SimpleAMGRuntime::SimpleAMGRuntime()
    : IterativeSolverRuntime(), mSetup(), mCurrentLevel( 0 ), mLibHandle( 0 ), mHostOnlyLevel(
          std::numeric_limits<IndexType>::max() ), mHostOnlyVars( 0 ), mReplicatedLevel(
          std::numeric_limits<IndexType>::max() )
{
}

SimpleAMG::~SimpleAMG()
{
}

SimpleAMG::SimpleAMGRuntime::~SimpleAMGRuntime()
{
}

void SimpleAMG::loadSetupLibs()
{
    std::string amgSetupLibrary;

    if ( common::Settings::getEnvironment( amgSetupLibrary, "SCAI_AMG_SETUP_LIBRARY" ) )
    {
        SCAI_LOG_INFO( logger, "Load all module libraries in " << amgSetupLibrary  )
        scai::common::LibModule::loadLibsByPath( amgSetupLibrary.c_str() );
    }
    else
    {
        SCAI_LOG_WARN( logger, "SCAI_AMG_SETUP_LIBRARY not set, take SingleGridSetup" )
    }
}

void SimpleAMG::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.SimpleAMG.initialize" )
    SCAI_LOG_DEBUG( logger, "initialize AMG, coefficients matrix = " << coefficients )
    SimpleAMGRuntime& runtime = getRuntime();

    if ( runtime.mSetup.get() == NULL )
    {
        loadSetupLibs();
    }

    // Info about available AMGSetup
    std::vector<std::string> values;  // string is create type for the factory
    AMGSetup::getCreateValues( values );

    SCAI_LOG_INFO( logger, "Factory of AMGSetup: " << values.size() << " entries" )

    for ( size_t i = 0; i < values.size(); ++i )
    {
        SCAI_LOG_DEBUG( logger, "   Registered values[" << i << "] = " << values[i] )
    }

    if ( runtime.mSetup.get() == NULL )
    {
        // no setup defined yet, so we take on from the factory
        if ( AMGSetup::canCreate( "SAMGPSetup" ) )
        {
            runtime.mSetup.reset( AMGSetup::create( "SAMGPSetup" ) );
            SCAI_LOG_INFO( logger, "SimpleAMG: take SAMGPSetup as AMGSetup" )
        }
        else if ( AMGSetup::canCreate( "SimpleAMGSetup" ) )
        {
            runtime.mSetup.reset( AMGSetup::create( "SimpleAMGSetup" ) );
            SCAI_LOG_INFO( logger, "SimpleAMG: take SimpleAMGSetup as AMGSetup" )
        }
        else if ( AMGSetup::canCreate( "SingleGridSetup" ) )
        {
            runtime.mSetup.reset( AMGSetup::create( "SingleGridSetup" ) );
            SCAI_LOG_INFO( logger, "SimpleAMG: take SingleGridSetup as AMGSetup" )
        }
    }

    if ( !runtime.mSetup )
    {
        COMMON_THROWEXCEPTION( "No AMGSetup found" )
    }

    AMGSetup& amgSetup = *runtime.mSetup;
    amgSetup.setMaxLevels( mMaxLevels );
    amgSetup.setMinVarsCoarseLevel( mMinVarsCoarseLevel );
    amgSetup.setHostOnlyLevel( runtime.mHostOnlyLevel );
    amgSetup.setHostOnlyVars( runtime.mHostOnlyVars );
    amgSetup.setReplicatedLevel( runtime.mReplicatedLevel );
    amgSetup.setCoarseLevelSolver( mCoarseLevelSolver );
    amgSetup.setSmoother( mSmoother );
    logSetupSettings();
    amgSetup.initialize( coefficients );
    logSetupInfo();
    logSolverInfo();
    logSetupDetails();

    if ( mSmootherContext )
    {
        for ( IndexType level = 0; level < ( IndexType ) amgSetup.getNumLevels() - 1; ++level )
        {
            amgSetup.getSmoother( level ).setContextPtr( mSmootherContext );
        }
    }

    IterativeSolver::initialize( coefficients );
    totalSmootherTime = 0.0;
    totalTransferTime = 0.0;
    totalResidualTime = 0.0;
    totalIterations = 0;
}

void SimpleAMG::iterate()
{
    SCAI_REGION( "Solver.SimpleAMG.iterate" )
    cycle();
    totalIterations++;
}

double SimpleAMG::getAverageSmootherTime() const
{
    return ( totalSmootherTime / totalIterations );
}

double SimpleAMG::getAverageTransferTime() const
{
    return ( totalTransferTime / totalIterations );
}

double SimpleAMG::getAverageResidualTime() const
{
    return ( totalResidualTime / totalIterations );
}

void SimpleAMG::setMaxLevels( unsigned int levels )
{
    mMaxLevels = levels;
}

void SimpleAMG::setMinVarsCoarseLevel( unsigned int vars )
{
    mMinVarsCoarseLevel = vars;
}

const Matrix& SimpleAMG::getGalerkin( unsigned int level )
{
    return getRuntime().mSetup->getGalerkin( level );
}

const Matrix& SimpleAMG::getRestriction( unsigned int level )
{
    return getRuntime().mSetup->getRestriction( level );
}

const Matrix& SimpleAMG::getInterpolation( unsigned int level )
{
    return getRuntime().mSetup->getInterpolation( level );
}

Vector& SimpleAMG::getSolutionVector( unsigned int level )
{
    return getRuntime().mSetup->getSolutionVector( level );
}

Vector& SimpleAMG::getRhsVector( unsigned int level )
{
    return getRuntime().mSetup->getRhsVector( level );
}

Solver& SimpleAMG::getSmoother( unsigned int level )
{
    return getRuntime().mSetup->getSmoother( level );
}

Solver& SimpleAMG::getCoarseLevelSolver()
{
    return getRuntime().mSetup->getCoarseLevelSolver();
}

void SimpleAMG::setSmootherContext( ContextPtr smootherContext )
{
    mSmootherContext = smootherContext;
}

void SimpleAMG::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    getRuntime().mHostOnlyLevel = hostOnlyLevel;
}

void SimpleAMG::setHostOnlyVars( IndexType hostOnlyVars )
{
    getRuntime().mHostOnlyVars = hostOnlyVars;
}

void SimpleAMG::setReplicatedLevel( IndexType replicatedLevel )
{
    getRuntime().mReplicatedLevel = replicatedLevel;
}

void SimpleAMG::setCoarseLevelSolver( SolverPtr solver )
{
    SCAI_LOG_DEBUG ( logger, "Set Coarse Level Solver to" << *solver )
    mCoarseLevelSolver = solver;
}

void SimpleAMG::setSmoother( SolverPtr solver )
{
    SCAI_LOG_DEBUG( logger, "Defined smoother for all level " << *solver )
    mSmoother = solver;
}

unsigned int SimpleAMG::getNumLevels()
{
    return getRuntime().mSetup->getNumLevels();
}

void SimpleAMG::cycle()
{
    // go via pointers because of const rhs on finest grid
    SimpleAMGRuntime& runtime = getRuntime();
    SCAI_REGION_N( "Solver.SimpleAMG.cycle", runtime.mCurrentLevel )
    // dereferences to current level solution + rhs
    common::shared_ptr<AMGSetup>& amgSetup = runtime.mSetup;
    const Vector* curRhsPtr = runtime.mRhs;
    Vector* curSolutionPtr = 0;

    if ( runtime.mCurrentLevel == 0 )
    {
        curSolutionPtr = &( runtime.mSolution.getReference() );
    }
    else
    {
        curSolutionPtr = &( amgSetup->getSolutionVector( runtime.mCurrentLevel ) );
        curRhsPtr = &( amgSetup->getRhsVector( runtime.mCurrentLevel ) );
    }

    Vector& curSolution = ( *curSolutionPtr );
    const Vector& curRhs = ( *curRhsPtr );

    //no more Smoothers we are on the coareste level
    if ( runtime.mCurrentLevel >= amgSetup->getNumLevels() - 1 )
    {
        amgSetup->getCoarseLevelSolver().solve( curSolution, curRhs );
    }
    else
    {
        const Matrix& curGalerkin = amgSetup->getGalerkin( runtime.mCurrentLevel );
        const Matrix& curRestriction = amgSetup->getRestriction( runtime.mCurrentLevel );
        const Matrix& curInterpolation = amgSetup->getInterpolation( runtime.mCurrentLevel );
        Vector& curTmpRhs = amgSetup->getTmpResVector( runtime.mCurrentLevel );
        Vector& curCoarseSolution = amgSetup->getSolutionVector( runtime.mCurrentLevel + 1 );
        Vector& curCoarseRhs = amgSetup->getRhsVector( runtime.mCurrentLevel + 1 );
        Solver& curSmoother = amgSetup->getSmoother( runtime.mCurrentLevel );
        // PreSmoothing
        SCAI_LOG_DEBUG( logger, "Pre smoothing on level " << runtime.mCurrentLevel )
        double smootherStartTime = common::Walltime::get();
        curSmoother.solve( curSolution, curRhs );
        totalSmootherTime += common::Walltime::get() - smootherStartTime;
        // Restrict residual to next coarser grid
        // and initialize solution
        SCAI_LOG_DEBUG( logger, "curTmpRhs=curRhs - curGalerkin * curSolution on level " << runtime.mCurrentLevel )
        double residualStartTime = common::Walltime::get();
        curTmpRhs = curRhs - curGalerkin * curSolution;
        totalResidualTime += common::Walltime::get() - residualStartTime;
        SCAI_LOG_DEBUG( logger, "curCoarseRhs = curRestriction * curTmpRhs on level " << runtime.mCurrentLevel )
        double transferStartTime = common::Walltime::get();
        curCoarseRhs = curRestriction * curTmpRhs;
        curCoarseSolution = 0.0;
        totalTransferTime += common::Walltime::get() - transferStartTime;
        ++runtime.mCurrentLevel;
        cycle();
        --runtime.mCurrentLevel;
        SCAI_LOG_DEBUG( logger,
                        "curSolution = curSolution + curInterpolation * curCoarseSolution on level " << runtime.mCurrentLevel )
        transferStartTime = common::Walltime::get();
        curSolution = curSolution + curInterpolation * curCoarseSolution;
        totalTransferTime += common::Walltime::get() - transferStartTime;
        SCAI_LOG_DEBUG( logger, "Post smoothing on level " << runtime.mCurrentLevel )
        smootherStartTime = common::Walltime::get();
        curSmoother.solve( curSolution, curRhs );
        totalSmootherTime += common::Walltime::get() - smootherStartTime;
    }
}

void SimpleAMG::logSetupSettings()
{
    if ( mLogger->getLogLevel() < LogLevel::solverInformation )
    {
        return;
    }

    if ( getRuntime().mSetup.get() != 0 )
    {
        mLogger->logMessage( LogLevel::solverInformation, "Running SimpleAMG.\n" );
        mLogger->logNewLine( LogLevel::solverInformation );
        mLogger->logMessage( LogLevel::solverInformation, "Setup Components:\n" );
        mLogger->logMessage( LogLevel::solverInformation, "=================\n" );
        std::stringstream outputCouplingPredicate;
        outputCouplingPredicate << getRuntime().mSetup->getCouplingPredicateInfo() << "\n";
        mLogger->logMessage( LogLevel::solverInformation, outputCouplingPredicate.str() );
        std::stringstream outputColoring;
        outputColoring << getRuntime().mSetup->getColoringInfo() << "\n";
        mLogger->logMessage( LogLevel::solverInformation, outputColoring.str() );
        std::stringstream outputInterpolation;
        outputInterpolation << getRuntime().mSetup->getInterpolationInfo() << "\n";
        mLogger->logMessage( LogLevel::solverInformation, outputInterpolation.str() );
        mLogger->logNewLine( LogLevel::solverInformation );
        std::stringstream outputRunning;
        outputRunning << "Running Setup (maxLevels=" << mMaxLevels << ", mMinVarsCoarseLevel=" << mMinVarsCoarseLevel
                      << ")...\n";
        mLogger->logMessage( LogLevel::solverInformation, outputRunning.str() );
    }
    else
    {
        mLogger->logMessage( LogLevel::solverInformation, "Running Gauss-Seidel.\n" );
    }
}

void SimpleAMG::logSetupInfo()
{
    if ( mLogger->getLogLevel() < LogLevel::solverInformation )
    {
        return;
    }

    mLogger->logMessage( LogLevel::solverInformation, "Setup done.\n" );
    mLogger->logType( LogLevel::solverInformation, "Number of levels\t: ", getRuntime().mSetup->getNumLevels() );
    mLogger->logNewLine( LogLevel::solverInformation );

    if ( mLogger->getLogLevel() < LogLevel::advancedInformation )
    {
        return;
    }

    for ( unsigned int i = 0; i < getRuntime().mSetup->getNumLevels(); ++i )
    {
        if ( i == 0 )
        {
            std::stringstream output1;
            output1 << "Operator Matrix Hierarchy:\n";
            mLogger->logMessage( LogLevel::advancedInformation, output1.str() );
            std::stringstream output2;
            output2 << "Lvl    #Rows    #Cols  #Entries Average\n";
            mLogger->logMessage( LogLevel::advancedInformation, output2.str() );
        }

        std::stringstream output;
        double averageNumValues = static_cast<double>( getRuntime().mSetup->getGalerkin( i ).getNumValues() )
                                  / getRuntime().mSetup->getGalerkin( i ).getNumRows();
        output << std::setw( 3 ) << i << std::setw( 9 ) << getRuntime().mSetup->getGalerkin( i ).getNumRows()
               << std::setw( 9 ) << getRuntime().mSetup->getGalerkin( i ).getNumColumns() << std::setw( 10 )
               << getRuntime().mSetup->getGalerkin( i ).getNumValues() << "  " << std::setw( 6 ) << std::fixed
               << std::setprecision( 1 ) << averageNumValues << "\n";
        mLogger->logMessage( LogLevel::advancedInformation, output.str() );
    }

    mLogger->logNewLine( LogLevel::advancedInformation );

    for ( unsigned int i = 0; i < getRuntime().mSetup->getNumLevels() - 1; ++i )
    {
        if ( i == 0 )
        {
            std::stringstream output1;
            output1 << "Interpolation Matrix Hierarchy:\n";
            mLogger->logMessage( LogLevel::advancedInformation, output1.str() );
            std::stringstream output2;
            output2 << "Lvl    #Rows    #Cols  #Entries Average\n";
            mLogger->logMessage( LogLevel::advancedInformation, output2.str() );
        }

        std::stringstream output;
        double averageNumValues = static_cast<double>( getRuntime().mSetup->getInterpolation( i ).getNumValues() )
                                  / getRuntime().mSetup->getInterpolation( i ).getNumRows();
        output << std::setw( 3 ) << i << std::setw( 9 ) << getRuntime().mSetup->getInterpolation( i ).getNumRows()
               << std::setw( 9 ) << getRuntime().mSetup->getInterpolation( i ).getNumColumns()
               << std::setw( 10 ) << getRuntime().mSetup->getInterpolation( i ).getNumValues() << "  "
               << std::setw( 6 ) << std::fixed << std::setprecision( 1 ) << averageNumValues << "\n";
        mLogger->logMessage( LogLevel::advancedInformation, output.str() );
    }

    mLogger->logNewLine( LogLevel::advancedInformation );
}

void SimpleAMG::logSolverInfo()
{
    if ( mLogger->getLogLevel() < LogLevel::solverInformation )
    {
        return;
    }

    if ( getRuntime().mSetup.get() != 0 )
    {
        mLogger->logNewLine( LogLevel::solverInformation );
        mLogger->logMessage( LogLevel::solverInformation, "Solution Components:\n" );
        mLogger->logMessage( LogLevel::solverInformation, "====================\n" );
        mLogger->logMessage( LogLevel::solverInformation, "Iteration Type : V-Cycle\n" );
        std::stringstream smootherInfo;
        smootherInfo << "Smoothing      : " << getRuntime().mSetup->getSmootherInfo() << "\n";
        mLogger->logMessage( LogLevel::solverInformation, smootherInfo.str() );
        std::stringstream coarseLevelSolverInfo;
        coarseLevelSolverInfo << "Coarse Level   : " << getRuntime().mSetup->getCoarseLevelSolverInfo() << "\n";
        mLogger->logMessage( LogLevel::solverInformation, coarseLevelSolverInfo.str() );
        mLogger->logNewLine( LogLevel::solverInformation );
    }
}

void SimpleAMG::logSetupDetails()
{
    if ( mLogger->getLogLevel() < LogLevel::advancedInformation )
    {
        return;
    }

    double sizeVector = 0.0;
    double sizeInterpolation = 0.0;
    double sizeRestriction = 0.0;
    double sizeGalerkin = 0.0;
    double sizeInterpolationCSR = 0.0;
    double sizeRestrictionCSR = 0.0;
    double sizeGalerkinCSR = 0.0;

    for ( unsigned int i = 0; i < getRuntime().mSetup->getNumLevels(); ++i )
    {
        // Vector
        if ( i == 0 )
        {
            sizeVector += 2 * getRuntime().mSetup->getGalerkin( 0 ).getValueTypeSize()
                          * getRuntime().mSetup->getGalerkin( 0 ).getNumRows();
        }
        else
        {
            sizeVector += getRuntime().mSetup->getSolutionVector( i ).getMemoryUsage();
            sizeVector += getRuntime().mSetup->getRhsVector( i ).getMemoryUsage();
        }

        if ( i != getRuntime().mSetup->getNumLevels() - 1 )
        {
            sizeVector += getRuntime().mSetup->getTmpResVector( i ).getMemoryUsage();
            // Interpolation
            {
                sizeInterpolation += getRuntime().mSetup->getInterpolation( i ).getMemoryUsage();
                IndexType numIndexInterpolationCSR = getRuntime().mSetup->getInterpolation( i ).getNumValues()
                                                     + getRuntime().mSetup->getInterpolation( i ).getNumRows();
                sizeInterpolationCSR += numIndexInterpolationCSR * sizeof( IndexType );
                IndexType numValueInterpolationCSR = getRuntime().mSetup->getInterpolation( i ).getNumValues();
                size_t interpolationSizeType = getRuntime().mSetup->getInterpolation( i ).getValueTypeSize();
                sizeInterpolationCSR += numValueInterpolationCSR * interpolationSizeType;
            }
            // Restriction
            {
                sizeRestriction += getRuntime().mSetup->getRestriction( i ).getMemoryUsage();
                IndexType numIndexRestrictionCSR = getRuntime().mSetup->getRestriction( i ).getNumValues()
                                                   + getRuntime().mSetup->getRestriction( i ).getNumRows();
                sizeRestrictionCSR += numIndexRestrictionCSR * sizeof( IndexType );
                size_t restrictionSizeType = getRuntime().mSetup->getRestriction( i ).getValueTypeSize();
                IndexType numValueRestrictionCSR = getRuntime().mSetup->getRestriction( i ).getNumValues();
                sizeRestrictionCSR += numValueRestrictionCSR * restrictionSizeType;
            }
        }

        // Galerkin
        {
            sizeGalerkin += getRuntime().mSetup->getGalerkin( i ).getMemoryUsage();
            IndexType numIndexGalerkinCSR = getRuntime().mSetup->getGalerkin( i ).getNumValues()
                                            + getRuntime().mSetup->getGalerkin( i ).getNumRows();
            sizeGalerkinCSR += numIndexGalerkinCSR * sizeof( IndexType );
            size_t galerkinSizeType = getRuntime().mSetup->getGalerkin( i ).getValueTypeSize();
            IndexType numValueGalerkinCSR = getRuntime().mSetup->getGalerkin( i ).getNumValues();
            sizeGalerkinCSR += numValueGalerkinCSR * galerkinSizeType;
        }
    }

    int overheadInterpolation = static_cast<int>( 100 * sizeInterpolation / sizeInterpolationCSR ) - 100;
    int overheadRestriction = static_cast<int>( 100 * sizeRestriction / sizeRestrictionCSR ) - 100;
    int overheadGalerkin = static_cast<int>( 100 * sizeGalerkin / sizeGalerkinCSR ) - 100;
    int overhead = static_cast<int>( 100 * ( sizeInterpolation + sizeRestriction + sizeGalerkin )
                                     / ( sizeInterpolationCSR + sizeRestrictionCSR + sizeGalerkinCSR ) ) - 100;
    size_t cgSolverValueTypeSize = getRuntime().mSetup->getGalerkin( getRuntime().mSetup->getNumLevels() - 1 ).getValueTypeSize();
    double sizeCGSolver = static_cast<double>( cgSolverValueTypeSize
                          * getRuntime().mSetup->getGalerkin( getRuntime().mSetup->getNumLevels() - 1 ).getNumRows()
                          * getRuntime().mSetup->getGalerkin( getRuntime().mSetup->getNumLevels() - 1 ).getNumRows() );
    double sizeTotal = static_cast<double>( sizeVector + sizeInterpolation + sizeRestriction + sizeGalerkin
                                            + sizeCGSolver );
    double sizeTotalCSR = static_cast<double>( sizeVector + sizeInterpolationCSR + sizeRestrictionCSR + sizeGalerkinCSR
                          + sizeCGSolver );
    std::stringstream outputVector;
    outputVector << "Vector         " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                 << sizeVector / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                 << sizeVector / ( 1024 * 1024 ) << " MB" << std::endl;
    std::stringstream outputInterpolation;
    outputInterpolation << "Interpolation  " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                        << sizeInterpolationCSR / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed
                        << std::setprecision( 1 ) << sizeInterpolation / ( 1024 * 1024 ) << " MB" << std::setw( 4 )
                        << overheadInterpolation << "% " << std::endl;
    std::stringstream outputRestriction;
    outputRestriction << "Restriction    " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                      << sizeRestrictionCSR / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed
                      << std::setprecision( 1 ) << sizeRestriction / ( 1024 * 1024 ) << " MB" << std::setw( 4 )
                      << overheadRestriction << "% " << std::endl;
    std::stringstream outputGalerkin;
    outputGalerkin << "Galerkin       " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                   << sizeGalerkinCSR / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed
                   << std::setprecision( 1 ) << sizeGalerkin / ( 1024 * 1024 ) << " MB" << std::setw( 4 )
                   << overheadGalerkin << "% " << std::endl;
    std::stringstream outputCGSolver;
    outputCGSolver << "Coarse Inverse " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                   << sizeCGSolver / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                   << sizeCGSolver / ( 1024 * 1024 ) << " MB" << std::endl;
    std::stringstream outputTotal;
    outputTotal << "Total          " << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                << sizeTotalCSR / ( 1024 * 1024 ) << " MB" << std::setw( 8 ) << std::fixed << std::setprecision( 1 )
                << sizeTotal / ( 1024 * 1024 ) << " MB" << std::setw( 4 ) << overhead << "% " << std::endl;
    mLogger->logNewLine( LogLevel::advancedInformation );
    mLogger->logMessage( LogLevel::advancedInformation, "Memory needed for AMG hierarchy:\n" );
    mLogger->logMessage( LogLevel::advancedInformation, "==========================================\n" );
    mLogger->logMessage( LogLevel::advancedInformation, "                    CSR(ref.)  Actual  pad\n" );
    mLogger->logMessage( LogLevel::advancedInformation, outputVector.str() );
    mLogger->logMessage( LogLevel::advancedInformation, outputInterpolation.str() );
    mLogger->logMessage( LogLevel::advancedInformation, outputRestriction.str() );
    mLogger->logMessage( LogLevel::advancedInformation, outputGalerkin.str() );
    mLogger->logMessage( LogLevel::advancedInformation, outputCGSolver.str() );
    mLogger->logMessage( LogLevel::advancedInformation, "------------------------------------------\n" );
    mLogger->logMessage( LogLevel::advancedInformation, outputTotal.str() );
    mLogger->logNewLine( LogLevel::advancedInformation );
}

SimpleAMG::SimpleAMGRuntime& SimpleAMG::getRuntime()
{
    return mSimpleAMGRuntime;
}

const SimpleAMG::SimpleAMGRuntime& SimpleAMG::getConstRuntime() const
{
    return mSimpleAMGRuntime;
}

SolverPtr SimpleAMG::copy()
{
    return SolverPtr( new SimpleAMG( *this ) );
}

void SimpleAMG::writeAt( std::ostream& stream ) const
{
    stream << "SimpleAMG ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string SimpleAMG::createValue()
{
    return "SimpleAMG";
}

Solver* SimpleAMG::create( const std::string name )
{
    return new SimpleAMG( name );
}

} /* end namespace solver */

} /* end namespace scai */
