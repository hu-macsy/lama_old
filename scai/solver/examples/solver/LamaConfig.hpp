/**
 * @file LamaConfig.hpp
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
#include <vector>
#include <locale>

/** Class that handles commonly used configuration values for LAMA
 *
 *  - Communicator for distributions
 *  - Context ( Host, GPU ) for matrices, vectors
 *  - Sparse matrix format ( csr, ell, jds, dia, coo ) + creator for it
 *  - Value type ( sp or float, dp or double )
 *  - number of threads used for CPU computing
 *  - block size: number of threads / block used for GPU computing
 *
 *  Very important: DO NOT USE THIS CLASS for a GLOBAL variable
 *
 *  ( so its constructor might be called before LAMA factories are available )
 *
 *  \code
 *      LamaConfig lamaconf;    // WRONG: constructur might fail
 *
 *      int main ()
 *      {
 *          LamaConfig lamaconf();   // RIGHT: lama has registered
 *      }
 *  \endcode
 */

class LamaConfig : public scai::common::Printable
{

public:

    /** Constructor with default settings. */

    LamaConfig();

    /** Destructor also might free communicator, context */

    ~LamaConfig();

    /** Getter for the specified matrix format, might default */

    scai::lama::Matrix::MatrixStorageFormat getFormat( ) const;

    /** Getter for the value type to be used */

    scai::common::scalar::ScalarType getValueType() const;

    /** Getter for the solver id */

    const char* getSolverName( ) const;

    /** get a new matrix of the specified matrix format and value type. */

    scai::lama::Matrix* getMatrix();

    scai::hmemo::ContextPtr getContextPtr() const
    {
        return mContext;
    }

    const scai::hmemo::Context& getContext() const
    {
        return *getContextPtr();
    }

    scai::dmemo::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const scai::dmemo::Communicator& getCommunicator() const
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
        return getMaxIter() != nIndex;
    }

    /** Get the maximal number of iterations. */

    IndexType getMaxIter() const;

    scai::solver::LogLevel::LogLevel getLogLevel() const;

    scai::lama::Matrix::SyncKind getCommunicationKind() const
    {
        return mCommunicationKind;
    }

    bool useMetis() const
    {
        return mUseMetis;
    }

    float getWeight() const;

    static inline void printHelp( const char* progName );

private:

    std::string mSolverName;   // name of solver, used for factory

    scai::lama::Matrix::MatrixStorageFormat mMatrixFormat;

    scai::hmemo::ContextPtr mContext;

    scai::lama::Matrix::SyncKind     mCommunicationKind;

    scai::common::scalar::ScalarType   mValueType;          // value type to use

    scai::dmemo::CommunicatorPtr      mComm;

    mutable IndexType mMaxIter;

    scai::solver::LogLevel::LogLevel   mLogLevel;

    bool                       mUseMetis;

    float mWeight;
};

/* ---------------------------------------------------------------------------- */

LamaConfig::LamaConfig()
{
    using scai::common::Settings;

    mCommunicationKind = scai::lama::Matrix::SYNCHRONOUS;
    mComm              = scai::dmemo::Communicator::getCommunicatorPtr();
    mMaxIter           = nIndex;
    mValueType         = scai::common::scalar::UNKNOWN;
    mUseMetis          = false;
    mWeight            = -1.0f;    // stands for undefined

    mContext           = scai::hmemo::Context::getContextPtr( );

    bool isSet;

    if ( Settings::getEnvironment( isSet, "SCAI_ASYNCHRONOUS" ) )
    {
        if ( isSet )
        {
            mCommunicationKind = scai::lama::Matrix::ASYNCHRONOUS;
        }
    }

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

    mWeight = 1.0;   // as default

    if ( scai::common::Settings::getEnvironment( val, "SCAI_WEIGHT" ) )
    {
        float weight;
        int narg = sscanf( val.c_str(), "%f", &weight );

        if ( narg > 0 && weight >= 0.0f )
        {
            mWeight = weight;
        }
        else
        {
            std::cout << "SCAI_WEIGHT=" << val << " illegal" <<  std::endl;
        }
    }

    mValueType = scai::common::TypeTraits<RealType>::stype;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        scai::common::scalar::ScalarType type = scai::common::str2ScalarType( val.c_str() );


        if ( type == scai::common::scalar::UNKNOWN )
        {
            std::cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << std::endl;
        }
        else if ( scai::lama::Vector::canCreate( scai::lama::VectorCreateKeyType( scai::lama::Vector::DENSE, type ) ) )
        {
            mValueType = type;
        }
        else
        {
            std::cout << "SCAI_TYPE=" << val << " known, but not supported for matrix/vector" << std::endl;
        }
    }

    if ( mContext->getType() == scai::hmemo::Context::CUDA )
    {
        mMatrixFormat = scai::lama::Matrix::ELL;
    }
    else
    {
        mMatrixFormat = scai::lama::Matrix::CSR;
    }

    if ( scai::common::Settings::getEnvironment( val, "SCAI_FORMAT" ) )
    {
        // check if we can create a matrix of this type

        scai::lama::Format::MatrixStorageFormat format = scai::lama::str2Format( val.c_str() );

        if ( format != scai::lama::Format::UNDEFINED )
        {
            mMatrixFormat = format;
        }
    }

    mSolverName = "CG";

    if ( scai::common::Settings::getEnvironment( val, "SCAI_SOLVER" ) )
    {
        // check if solver is available

        if ( !scai::solver::Solver::canCreate( val ) )
        {
            std::cout << "ATTENTION: solver " << val << " not available" << std::endl;
        }
        else
        {
            mSolverName = val;
        }
    }

    // solver log level, default is convergence History

    mLogLevel = scai::solver::LogLevel::convergenceHistory;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_SOLVER_LOG" ) )
    {
        scai::solver::LogLevel::LogLevel level = scai::solver::str2LogLevel( val.c_str() );

        if ( level == scai::solver::LogLevel::UNKNOWN )
        {
            std::cout << "SCAI_SOLVER_LOG: " << val << " is not a logging level" << std::endl;
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
    using scai::common::Settings;

    stream << "LAMA configuration" << std::endl;
    stream << "==================" << std::endl;
    stream << "Solver            = " << mSolverName << std::endl;
    stream << "Solver Logging    = " << mLogLevel << std::endl;
    stream << "Context           = " << getContext() << std::endl;
    stream << "Communicator      = " << *mComm << std::endl;
    stream << "Matrix format     = " << getFormat() << std::endl;
    stream << "CommKind          = " << mCommunicationKind << std::endl;
    stream << "ValueType         = " << getValueType() << std::endl;
    stream << "#Threads/CPU      = " << omp_get_max_threads() << std::endl;
    stream << "weight            = " << getWeight() << std::endl;

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
}

scai::lama::Matrix::MatrixStorageFormat LamaConfig::getFormat( ) const
{
    return mMatrixFormat;
}

scai::common::scalar::ScalarType LamaConfig::getValueType() const
{
    return mValueType;
}

scai::solver::LogLevel::LogLevel LamaConfig::getLogLevel() const
{
    return mLogLevel;
}

const char* LamaConfig::getSolverName( ) const
{
    return mSolverName.c_str();
}

float LamaConfig::getWeight() const
{
    return mWeight;
}

IndexType LamaConfig::getMaxIter() const
{
    if ( mMaxIter == nIndex )
    {
        // not defined yet

        IndexType iter;

        if ( scai::common::Settings::getEnvironment( iter, "SCAI_MAX_ITER" ) )
        {
            mMaxIter = iter;
        }
    }

    return mMaxIter;
}

scai::lama::Matrix* LamaConfig::getMatrix()
{
    return scai::lama::Matrix::getMatrix( mMatrixFormat, mValueType );
}

void LamaConfig::printHelp( const char* progName )
{
    using std::cout;
    using std::endl;

    cout << "Usage: " << progName << " <filename> [ options ][ T<num_threads> ] " << endl;
    cout << "         solver specific options:" << endl;
    cout << "         --SCAI_SOLVER=[CG|BiCG|...]" << endl;
    cout << "         --SCAI_SOLVER_LOG=[noLogging|convergenceHistory|solverInformation|advancedInformation|completeInformation]" << endl;
    cout << "         --SCAI_MAX_ITER=<int_val>" << endl;
    cout << "         --SCAI_FORMAT=[CSR|ELL|JDS|DIA|COO]" << endl;
    cout << "         --SCAI_TYPE=[float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLongDouble]" << endl;
    cout << "         --SCAI_NUM_THREADS=..." << endl;
    cout << "         --SCAI_USE_METIS=<flag>" << endl;
    cout << "         --SCAI_ASYNCHRONOUS=<flag>" << endl;
    cout << "         or general options:" << endl;
    cout << "         --SCAI_COMMUNICATOR=[MPI|GPI|NO]" << endl;
    cout << "         --SCAI_CONTEXT=[Host|CUDA|MIC]" << endl;
    cout << "         --SCAI_DEVICE=[0|1|...]" << endl;
    cout << "         --SCAI_CUDA_USE_TEXTURE=[0|1]" << endl;
    cout << "         --SCAI_CUDA_USE_SHARED_MEM=[0|1]" << endl;
    cout << "         --SCAI_CUDA_BLOCK_SIZE=[64|128|...]" << endl;
}

