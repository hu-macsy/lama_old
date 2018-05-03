/**
 * @file LamaConfig.hpp
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
 * @brief Structure that contains configuration for LAMA
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include <scai/lama.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/common/Printable.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/solver/Solver.hpp>
#include <scai/solver/logger/LogLevel.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/Settings.hpp>

#include <cstring>
#include <sstream>
#include <vector>
#include <locale>

namespace scai
{

/** Class that handles commonly used configuration values for LAMA
 *
 *  - Communicator for distributions
 *  - Context ( Host, GPU ) for matrices, vectors
 *  - Sparse matrix format ( csr, ell, jds, dia, coo ) + creator for it
 *  - Value type ( float, double, ... )
 *  - number of threads used for CPU computing
 *  - block size: number of threads / block used for GPU computing
 *
 *  Very important: DO NOT USE THIS CLASS for a GLOBAL variable
 *  ( so its constructor might be called before LAMA factories are available )
 *
 *  Very important: Call constructor after command line arguments have been parsed
 *  ( otherwise only envrionment variables will take effect )
 *
 *  \code
 *      LamaConfig lamaconf;    // WRONG: constructur might fail
 *
 *      int main( int argc, const char* argv[] )
 *      {
 *           common::Settings::parseArgs( argc, argv );
 *           LamaConfig lamaconf;                        // must be defined after parseArgs
 *           ....
 *      }
 *  \endcode
 */

class LamaConfig : public common::Printable
{

public:

    /** Constructor with default settings. */

    LamaConfig();

    /** Destructor also might free communicator, context */

    ~LamaConfig();

    /** Getter for the specified matrix format, might default */

    lama::Format getFormat( ) const;

    /** Getter for the value type to be used */

    common::ScalarType getValueType() const;

    /** Getter for the solver id */

    const std::string& getSolverName( ) const;

    /** Getter for the norm id */

    const std::string& getNorm( ) const;

    /** get a new matrix of the specified matrix format and value type. */

    template<typename ValueType>
    lama::Matrix<ValueType>* getMatrix();

    hmemo::ContextPtr getContextPtr() const
    {
        return mContext;
    }

    const hmemo::Context& getContext() const
    {
        return *getContextPtr();
    }

    dmemo::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const dmemo::Communicator& getCommunicator() const
    {
        return *mComm;
    }

    /** Implements writeAt for class Printable so you can use it as follows:
     *
     *  \code
     *  LamaConfig lamaconfig;
     *  ...
     *  std::cout << "Config = " << lamaconfig << std::endl;
     *  \endcode
     */

    void writeAt( std::ostream& stream ) const;

    /** Query if use has set maximal number of iterations. */

    bool hasMaxIter() const
    {
        return getMaxIter() != invalidIndex;
    }

    /** Get the maximal number of iterations. */

    IndexType getMaxIter() const;

    solver::LogLevel getLogLevel() const;

    lama::SyncKind getCommunicationKind() const
    {
        return mCommunicationKind;
    }

    bool useMetis() const
    {
        return mUseMetis;
    }

    double getAbsoluteTolerance() const
    {
        return mAbsoluteTolerance;
    }

    double getRelativeTolerance() const
    {
        return mRelativeTolerance;
    }

    double getDivergenceTolerance() const
    {
        return mDivergenceTolerance;
    }

    static inline void printHelp( const char* progName );

private:

    std::string mSolverName;   // name of solver, used for factory
    std::string mNorm;         // name of norm, not yet factory

    lama::Format mMatrixFormat;

    hmemo::ContextPtr mContext;

    lama::SyncKind    mCommunicationKind;

    common::ScalarType   mValueType;          // value type to use

    dmemo::CommunicatorPtr      mComm;

    IndexType mMaxIter;

    solver::LogLevel   mLogLevel;

    bool mUseMetis;

    double mRelativeTolerance;
    double mAbsoluteTolerance;
    double mDivergenceTolerance;
};

/* ---------------------------------------------------------------------------- */

#define CONFIG_ERROR( msg )                                    \
    {                                                              \
        std::ostringstream errorStr;                               \
        errorStr << msg;                                           \
        std::cout << "ERROR: " << errorStr.str()  <<  std::endl;   \
        LamaConfig::printHelp( "lamaSolver" );                     \
        exit( 1 );                                                 \
    }

/* ---------------------------------------------------------------------------- */

/** This help routine reads a non-negative double value form an environment variable */

static void getTolerance( double& tolerance, const char* name )
{
    tolerance = 0.0;  // as default

    std::string stringVal;

    if ( common::Settings::getEnvironment( stringVal, name ) )
    {
        std::istringstream input( stringVal );

        double tol;

        if ( input >> tol )
        {
            if ( tol >= 0.0 )
            {
                tolerance = tol;
            }
            else
            {
                CONFIG_ERROR( name << "="  << stringVal << " illegal, mus be positive" )
            }
        }
        else
        {
            CONFIG_ERROR( name << stringVal << " illegal, no double value" )
        }
    }
}

/* ---------------------------------------------------------------------------- */

LamaConfig::LamaConfig()
{
    using common::Settings;
 
    // take default communicator, can be set by SCAI_COMMUNICATOR

    mComm = dmemo::Communicator::getCommunicatorPtr();

    // allow settings specified for each process
    // Be careful: can cause problems, e.g. MAX_ITER, xxx_TOL

    common::Settings::setRank( mComm->getRank() );

    // take default context, can be set by SCAI_CONTEXT

    mContext = hmemo::Context::getContextPtr( );

    bool isSet;
    int kind;

    mCommunicationKind = lama::SyncKind::SYNCHRONOUS;

    if ( Settings::getEnvironment( kind, "SCAI_ASYNCHRONOUS" ) )
    {
        // 0: SYNC, 1 : ASYNC_COMM, 2 : ASYNC_LOCAL
        mCommunicationKind = lama::SyncKind( kind );
    }

    mUseMetis = false;  // default

    if ( Settings::getEnvironment( isSet, "SCAI_USE_METIS" ) )
    {
        mUseMetis = isSet;
    }

    int nThreads;

    if ( Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    std::string val;

    getTolerance( mRelativeTolerance, "SCAI_REL_TOL" );
    getTolerance( mAbsoluteTolerance, "SCAI_ABS_TOL" );
    getTolerance( mDivergenceTolerance, "SCAI_DIV_TOL" );

    mMaxIter = invalidIndex;

    common::Settings::getEnvironment( mMaxIter, "SCAI_MAX_ITER" );

    // ValueType to be used for vector/matrix

    mValueType = common::TypeTraits<DefaultReal>::stype;

    if ( common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        common::ScalarType type = common::str2ScalarType( val.c_str() );

        if ( type == common::ScalarType::UNKNOWN )
        {
            CONFIG_ERROR( "SCAI_TYPE=" << val << " illegal, is not a scalar type" )
        }
        else if ( lama::_Vector::canCreate( lama::VectorCreateKeyType( lama::VectorKind::DENSE, type ) ) )
        {
            mValueType = type;
        }
        else
        {
            CONFIG_ERROR( "SCAI_TYPE=" << val << " known, but not supported for matrix/vector" )
        }
    }

    if ( mContext->getType() == common::ContextType::CUDA )
    {
        mMatrixFormat = lama::Format::ELL;
    }
    else
    {
        mMatrixFormat = lama::Format::CSR;
    }

    if ( common::Settings::getEnvironment( val, "SCAI_FORMAT" ) )
    {
        // check if we can create a matrix of this type

        lama::Format format = lama::str2Format( val.c_str() );

        if ( format != lama::Format::UNDEFINED )
        {
            mMatrixFormat = format;
        }
    }

    mSolverName = "CG";

    if ( common::Settings::getEnvironment( val, "SCAI_SOLVER" ) )
    {
        // check if solver is available

        if ( !solver::_Solver::canCreate( solver::SolverCreateKeyType( mValueType, val ) ) )
        {
            CONFIG_ERROR( "solver " << val << " not available" )
        }
        else
        {
            mSolverName = val;
        }
    }

    mNorm = "L2";

    common::Settings::getEnvironment( mNorm, "SCAI_NORM" );

    if ( ! lama::Norm<DefaultReal>::canCreate( mNorm ) )
    {
        CONFIG_ERROR( "norm " << mNorm << " not available" )
    }

    // solver log level, default is convergence History

    mLogLevel = solver::LogLevel::convergenceHistory;

    if ( common::Settings::getEnvironment( val, "SCAI_SOLVER_LOG" ) )
    {
        solver::LogLevel level = solver::str2LogLevel( val.c_str() );

        if ( level == solver::LogLevel::UNKNOWN )
        {
            CONFIG_ERROR( "SCAI_SOLVER_LOG: " << val << " is not a logging level" )
        }
        else
        {
            mLogLevel = level;
        }
    }
}

LamaConfig::~LamaConfig()
{
}

void LamaConfig::writeAt( std::ostream& stream ) const
{
    using common::Settings;

    stream << "lamaSolver configuration" << std::endl;
    stream << "========================" << std::endl;
    stream << "Solver            = " << mSolverName << std::endl;
    stream << "Solver Logging    = " << mLogLevel << std::endl;
    stream << "Context           = " << getContext() << std::endl;
    stream << "Communicator      = " << *mComm << std::endl;
    stream << "Matrix format     = " << getFormat() << std::endl;
    stream << "CommKind          = " << mCommunicationKind << std::endl;
    stream << "ValueType         = " << getValueType() << std::endl;
    stream << "#Threads/CPU      = " << omp_get_max_threads() << std::endl;

    bool useTexture;
    bool useSharedMem;
    int  blockSize;

    if ( Settings::getEnvironment( useTexture, "SCAI_CUDA_USE_TEXTURE" ) )
    {
        stream << "useTexture(GPU)   = " << useTexture << std::endl;
    }

    if ( Settings::getEnvironment( useSharedMem, "SCAI_CUDA_USE_SHARED_MEM" ) )
    {
        stream << "useSharedMem(GPU) = " << useSharedMem << std::endl;
    }

    if ( Settings::getEnvironment( blockSize, "SCAI_CUDA_BLOCK_SIZE" ) )
    {
        stream << "BlockSize(GPU)    = " << blockSize << std::endl;
    }

    if ( hasMaxIter () )
    {
        stream << "Max iter          = " << mMaxIter << std::endl;
    }

    stream << "Tolerances        =";

    if ( mRelativeTolerance > 0.0 )
    {
        stream << " " << mRelativeTolerance << "(relative)";
    }

    if ( mAbsoluteTolerance > 0.0 )
    {
        stream << " " << mAbsoluteTolerance << "(absolute)";
    }

    if ( mDivergenceTolerance > 0.0 )
    {
        stream << " " << mDivergenceTolerance << "(divergence)";
    }

    stream << std::endl;

    stream << "Norm              = " << getNorm() << std::endl;
}

lama::Format LamaConfig::getFormat( ) const
{
    return mMatrixFormat;
}

common::ScalarType LamaConfig::getValueType() const
{
    return mValueType;
}

solver::LogLevel LamaConfig::getLogLevel() const
{
    return mLogLevel;
}

const std::string& LamaConfig::getSolverName( ) const
{
    return mSolverName;
}

const std::string& LamaConfig::getNorm( ) const
{
    return mNorm;
}

IndexType LamaConfig::getMaxIter() const
{
    return mMaxIter;
}

template<typename ValueType>
lama::Matrix<ValueType>* LamaConfig::getMatrix()
{
    return lama::Matrix<ValueType>::getMatrix( mMatrixFormat );
}

static std::string getLoggers()
{
    std::ostringstream loggerNames;

    std::vector<std::string> vals;
    solver::Solver<DefaultReal>::getCreateValues( vals );

    for ( size_t i = 0; i < vals.size(); ++i )
    {
        if ( i > 0 )
        {
            loggerNames << "|";
        }

        loggerNames << vals[i];
    }

    return loggerNames.str();
}

static std::string getLogLevels()
{
    std::ostringstream levelNames;

    for ( int i = 0; i < static_cast<int>( solver::LogLevel::UNKNOWN ); ++i )
    {
        if ( i > 0 )
        {
            levelNames << "|";
        }

        levelNames << solver::LogLevel( i );
    }

    return levelNames.str();
}

void LamaConfig::printHelp( const char* progName )
{
    using std::cout;
    using std::endl;

    cout << "Usage: " << progName << " <matrix_filename> [ rhs ] [start_solution] [ final_solution_filename ] [ options ] " << endl;
    cout << "         rhs, start_solution can be either a value or a filename" << endl;
    cout << "         solver specific options:" << endl;
    cout << "         --SCAI_DISTRIBUTION=<dist_filename>" << endl;
    cout << "         --SCAI_PARTITIONING=OFF|BLOCK|CYCLIC|METIS|PARAMETS" << endl;
    cout << "         --SCAI_SOLVER=" << getLoggers() << endl;
    cout << "         --SCAI_SOLVER_LOG=" << getLogLevels() << endl;
    cout << "         --SCAI_MAX_ITER=<int_val>" << endl;
    cout << "         --SCAI_NORM=L1|L2|Max" << endl;
    cout << "         --SCAI_REL_TOL=<val>" << endl;
    cout << "         --SCAI_ABS_TOL=<val>" << endl;
    cout << "         --SCAI_DIV_TOL=<val>" << endl;
    cout << "         --SCAI_FORMAT=[CSR|ELL|JDS|DIA|COO]" << endl;
    cout << "         --SCAI_TYPE=[float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLongDouble]" << endl;
    cout << "         --SCAI_NUM_THREADS=..." << endl;
    cout << "         --SCAI_USE_METIS=<flag>" << endl;
    cout << "         --SCAI_ASYNCHRONOUS=<flag>" << endl;
    cout << "         or general options:" << endl;
    cout << "         --SCAI_COMMUNICATOR=[MPI|NO]" << endl;
    cout << "         --SCAI_CONTEXT=[Host|CUDA]" << endl;
    cout << "         --SCAI_DEVICE=[0|1|...]" << endl;
    cout << "         --SCAI_CUDA_USE_TEXTURE=[0|1]" << endl;
    cout << "         --SCAI_CUDA_USE_SHARED_MEM=[0|1]" << endl;
    cout << "         --SCAI_CUDA_BLOCK_SIZE=[64|128|...]" << endl;
}

}
