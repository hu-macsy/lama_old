/**
 * @file SimpleAMG.cpp
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
#include <scai/common/macros/instantiate.hpp>

// std
#include <cstdlib>
#include <iomanip>

using namespace scai::hmemo;

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SimpleAMG<ValueType>::logger, 
                              "Solver.IterativeSolver.SimpleAMG" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* SimpleAMG<ValueType>::create()
{
    return new SimpleAMG<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType SimpleAMG<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "SimpleAMG" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
SimpleAMG<ValueType>::SimpleAMG( const std::string& id ) : 

    IterativeSolver<ValueType>( id ), 
    mMaxLevels( 25 ), 
    mMinVarsCoarseLevel( 100 )
{
    SCAI_LOG_INFO( logger, "SimpleAMG, id = " << id << " created, no logger" )
}

template<typename ValueType>
SimpleAMG<ValueType>::SimpleAMG( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id, logger ), 
    mMaxLevels( 25 ), 
    mMinVarsCoarseLevel( 100 )
{
    SCAI_LOG_INFO( SimpleAMG<ValueType>::logger, "SimpleAMG, id = " << id << " created, with logger" )
}

template<typename ValueType>
SimpleAMG<ValueType>::SimpleAMG( const SimpleAMG& other ) : 

    IterativeSolver<ValueType>( other ), 
    mMaxLevels( other.mMaxLevels ), 
    mMinVarsCoarseLevel( other.mMinVarsCoarseLevel )
{
}

template<typename ValueType>
SimpleAMG<ValueType>::SimpleAMGRuntime::SimpleAMGRuntime() : 

    IterativeSolver<ValueType>::IterativeSolverRuntime(), 
    mSetup(), 
    mCurrentLevel( 0 ), 
    mLibHandle( 0 ), 
    mHostOnlyLevel( std::numeric_limits<IndexType>::max() ), 
    mHostOnlyVars( 0 ), 
    mReplicatedLevel( std::numeric_limits<IndexType>::max() )
{
}

template<typename ValueType>
SimpleAMG<ValueType>::SimpleAMGRuntime::~SimpleAMGRuntime() 
{
    SCAI_LOG_INFO( logger, "~SimpleAMGRuntime, delete setup" )
    mSetup.reset();
    SCAI_LOG_INFO( logger, "~SimpleAMGRuntime" )
}

template<typename ValueType>
SimpleAMG<ValueType>::~SimpleAMG()
{
    SCAI_LOG_INFO( logger, "~SimpleAMG" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void SimpleAMG<ValueType>::loadSetupLibs()
{
    std::string amgSetupLibrary;

    if ( common::Settings::getEnvironment( amgSetupLibrary, "SCAI_LIBRARY_PATH" ) )
    {
        SCAI_LOG_INFO( logger, "Load all module libraries in " << amgSetupLibrary  )
        scai::common::LibModule::loadLibsByPath( amgSetupLibrary.c_str() );
    }
    else
    {
        SCAI_LOG_WARN( logger, "SCAI_LIBRARY_PATH not set, only defaults can be used" )
    }
}

template<typename ValueType>
void SimpleAMG<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.SimpleAMG.initialize" )
    SCAI_LOG_DEBUG( logger, "initialize AMG, coefficients matrix = " << coefficients )

    SimpleAMGRuntime& runtime = getRuntime();

    if ( runtime.mSetup.get() == NULL )
    {
        loadSetupLibs();
    }

    std::string amgSetupKey = "SingleGridSetup";    // take this as default

    if ( common::Settings::getEnvironment( amgSetupKey, "SCAI_AMG_SETUP" ) )
    {
        // give an error message if key is unknown, printing all possible keys

        if ( !AMGSetup<ValueType>::canCreate( amgSetupKey ) )
        {
            std::vector<std::string> values;

            AMGSetup<ValueType>::getCreateValues( values );

            std::string valuesStr;

            for ( auto const& v : values )
            {
                if ( valuesStr.length() )
                {
                    valuesStr += ":";
                }
                valuesStr += v;
            }

            SCAI_LOG_ERROR( logger, "SCAI_AMG_SETUP=" << amgSetupKey << ", key not available, only " << valuesStr )

            amgSetupKey = "SingleGridSetup";
        }
    }
    else
    {
        SCAI_LOG_WARN( logger, "Environment variable SCAI_AMG_SETUP not set, take default: " << amgSetupKey )
    }

    runtime.mSetup.reset( AMGSetup<ValueType>::getAMGSetup( amgSetupKey ) );

    if ( !runtime.mSetup )
    {
        COMMON_THROWEXCEPTION( "No AMGSetup found, key = " << amgSetupKey )
    }

    AMGSetup<ValueType>& amgSetup = *runtime.mSetup;
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

    IterativeSolver<ValueType>::initialize( coefficients );
    totalSmootherTime = 0.0;
    totalTransferTime = 0.0;
    totalResidualTime = 0.0;
    totalIterations = 0;
}

template<typename ValueType>
void SimpleAMG<ValueType>::iterate()
{
    SCAI_REGION( "Solver.SimpleAMG.iterate" )
    cycle();
    totalIterations++;
}

template<typename ValueType>
double SimpleAMG<ValueType>::getAverageSmootherTime() const
{
    return ( totalSmootherTime / totalIterations );
}

template<typename ValueType>
double SimpleAMG<ValueType>::getAverageTransferTime() const
{
    return ( totalTransferTime / totalIterations );
}

template<typename ValueType>
double SimpleAMG<ValueType>::getAverageResidualTime() const
{
    return ( totalResidualTime / totalIterations );
}

template<typename ValueType>
void SimpleAMG<ValueType>::setMaxLevels( IndexType levels )
{
    mMaxLevels = levels;
}

template<typename ValueType>
void SimpleAMG<ValueType>::setMinVarsCoarseLevel( IndexType vars )
{
    mMinVarsCoarseLevel = vars;
}

template<typename ValueType>
const Matrix<ValueType>& SimpleAMG<ValueType>::getGalerkin( IndexType level )
{
    return getRuntime().mSetup->getGalerkin( level );
}

template<typename ValueType>
const Matrix<ValueType>& SimpleAMG<ValueType>::getRestriction( IndexType level )
{
    return getRuntime().mSetup->getRestriction( level );
}

template<typename ValueType>
const Matrix<ValueType>& SimpleAMG<ValueType>::getInterpolation( IndexType level )
{
    return getRuntime().mSetup->getInterpolation( level );
}

template<typename ValueType>
Vector<ValueType>& SimpleAMG<ValueType>::getSolutionVector( IndexType level )
{
    return getRuntime().mSetup->getSolutionVector( level );
}

template<typename ValueType>
Vector<ValueType>& SimpleAMG<ValueType>::getRhsVector( IndexType level )
{
    return getRuntime().mSetup->getRhsVector( level );
}

template<typename ValueType>
Solver<ValueType>& SimpleAMG<ValueType>::getSmoother( IndexType level )
{
    return getRuntime().mSetup->getSmoother( level );
}

template<typename ValueType>
Solver<ValueType>& SimpleAMG<ValueType>::getCoarseLevelSolver()
{
    return getRuntime().mSetup->getCoarseLevelSolver();
}

template<typename ValueType>
void SimpleAMG<ValueType>::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    getRuntime().mHostOnlyLevel = hostOnlyLevel;
}

template<typename ValueType>
void SimpleAMG<ValueType>::setHostOnlyVars( IndexType hostOnlyVars )
{
    getRuntime().mHostOnlyVars = hostOnlyVars;
}

template<typename ValueType>
void SimpleAMG<ValueType>::setReplicatedLevel( IndexType replicatedLevel )
{
    getRuntime().mReplicatedLevel = replicatedLevel;
}

template<typename ValueType>
void SimpleAMG<ValueType>::setCoarseLevelSolver( SolverPtr<ValueType> solver )
{
    SCAI_LOG_DEBUG ( logger, "Set Coarse Level Solver to" << *solver )
    mCoarseLevelSolver = solver;
}

template<typename ValueType>
void SimpleAMG<ValueType>::setSmoother( SolverPtr<ValueType> solver )
{
    SCAI_LOG_DEBUG( logger, "Defined smoother for all level " << *solver )
    mSmoother = solver;
}

template<typename ValueType>
IndexType SimpleAMG<ValueType>::getNumLevels()
{
    return getRuntime().mSetup->getNumLevels();
}

template<typename ValueType>
void SimpleAMG<ValueType>::cycle()
{
    // go via pointers because of const rhs on finest grid
    SimpleAMGRuntime& runtime = getRuntime();
    SCAI_REGION_N( "Solver.SimpleAMG.cycle", runtime.mCurrentLevel )
    // dereferences to current level solution + rhs
    AMGSetup<ValueType>& amgSetup = *runtime.mSetup;
    const Vector<ValueType>* curRhsPtr = runtime.mRhs;
    Vector<ValueType>* curSolutionPtr = 0;

    if ( runtime.mCurrentLevel == 0 )
    {
        curSolutionPtr = &( runtime.mSolution.getReference() );
    }
    else
    {
        curSolutionPtr = &( amgSetup.getSolutionVector( runtime.mCurrentLevel ) );
        curRhsPtr = &( amgSetup.getRhsVector( runtime.mCurrentLevel ) );
    }

    Vector<ValueType>& curSolution = ( *curSolutionPtr );
    const Vector<ValueType>& curRhs = ( *curRhsPtr );

    //no more Smoothers we are on the coareste level
    if ( runtime.mCurrentLevel >= amgSetup.getNumLevels() - 1 )
    {
        amgSetup.getCoarseLevelSolver().solve( curSolution, curRhs );
    }
    else
    {
        const Matrix<ValueType>& curGalerkin = amgSetup.getGalerkin( runtime.mCurrentLevel );
        const Matrix<ValueType>& curRestriction = amgSetup.getRestriction( runtime.mCurrentLevel );
        const Matrix<ValueType>& curInterpolation = amgSetup.getInterpolation( runtime.mCurrentLevel );
        Vector<ValueType>& curTmpRhs = amgSetup.getTmpResVector( runtime.mCurrentLevel );
        Vector<ValueType>& curCoarseSolution = amgSetup.getSolutionVector( runtime.mCurrentLevel + 1 );
        Vector<ValueType>& curCoarseRhs = amgSetup.getRhsVector( runtime.mCurrentLevel + 1 );
        Solver<ValueType>& curSmoother = amgSetup.getSmoother( runtime.mCurrentLevel );
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
        curCoarseSolution.setSameValue( curInterpolation.getColDistributionPtr(), 0 );
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

template<typename ValueType>
void SimpleAMG<ValueType>::logSetupSettings()
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

template<typename ValueType>
void SimpleAMG<ValueType>::logSetupInfo()
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

    for ( IndexType i = 0; i < getRuntime().mSetup->getNumLevels(); ++i )
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

    for ( IndexType i = 0; i < getRuntime().mSetup->getNumLevels() - 1; ++i )
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

template<typename ValueType>
void SimpleAMG<ValueType>::logSolverInfo()
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

template<typename ValueType>
void SimpleAMG<ValueType>::logSetupDetails()
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

    for ( IndexType i = 0; i < getRuntime().mSetup->getNumLevels(); ++i )
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

template<typename ValueType>
typename SimpleAMG<ValueType>::SimpleAMGRuntime& SimpleAMG<ValueType>::getRuntime()
{
    return mSimpleAMGRuntime;
}

template<typename ValueType>
const typename SimpleAMG<ValueType>::SimpleAMGRuntime& SimpleAMG<ValueType>::getRuntime() const
{
    return mSimpleAMGRuntime;
}

template<typename ValueType>
SimpleAMG<ValueType>* SimpleAMG<ValueType>::copy()
{
    return new SimpleAMG( *this );
}

template<typename ValueType>
void SimpleAMG<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "SimpleAMG<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SimpleAMG, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
