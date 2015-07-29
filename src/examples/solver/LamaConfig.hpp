/**
 * @file LamaConfig.hpp
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
 * @brief Structure that contains configuration for LAMA
 * @author Thomas Brandes
 * @date 10.06.2013
 * @since 1.1.0
 */

#pragma once

#include <lama.hpp>

#include <memory/Context.hpp>
#include <common/Printable.hpp>
#include <lama/matrix/all.hpp>
#include <lama/solver/logger/LogLevel.hpp>
#include <omp.h>

#include <cstring>
#include <vector>

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

class LamaConfig : public Printable
{

public:

    /** Constructor with default settings. */

    LamaConfig();

    /** Destructor also might free communicator, context */

    ~LamaConfig();

    /** set an argument for configuration, e.g.  host, cuda, gpu, ell */

    void setArg( const char* arg );

    /** Getter for the specified matrix format, might default */

    const char* getFormat( ) const;

    /** Create a new sparse matrix of the desired matrix type. */

    template<typename ValueType>
    lama::SparseMatrix<ValueType>* createSparseMatrix();

    template<typename ValueType>
    lama::SparseMatrix<ValueType>* createSparseMatrix( const char* format );

    memory::ContextPtr getContextPtr() const
    {
        // Create a new context if not done yet

        if ( !mContext )
        {
            mContext = memory::Context::getContextPtr( mContextType, mDevice );
        }

        return mContext;
    }

    const memory::Context& getContext() const
    {
        return *getContextPtr();
    }

    lama::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const lama::Communicator& getCommunicator() const
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
        return mMaxIter != lama::nIndex; 
    }

    /** Get the maximal number of iterations. */

    lama::IndexType getMaxIter() const
    { 
        return mMaxIter; 
    }

    lama::LogLevel::LogLevel getLogLevel() const
    {
        return mLogLevel;
    }

    lama::Matrix::SyncKind getCommunicationKind() const
    {
        return mCommunicationKind;
    }

    common::ScalarType getValueType() const
    {
        return mValueType;
    }

    bool useMetis() const
    {
        return mUseMetis;
    }

    float getWeight() const
    {
        return mWeight;
    }

private:

    std::string              mMatrixFormat;

    memory::ContextType        mContextType;

    mutable memory::ContextPtr   mContext;

    lama::Matrix::SyncKind     mCommunicationKind;

    common::ScalarType   mValueType;          // value type to use

    lama::CommunicatorPtr      mComm;

    lama::IndexType            mMaxIter;

    lama::LogLevel::LogLevel   mLogLevel;

    bool                       mUseMetis;

    float                      mWeight;

    int                        mDevice;

    /** Help routine to query if argument has only digits. */

    inline bool isNumber( const char* arg );

    inline bool isReal( const char* arg );
};

/* ---------------------------------------------------------------------------- */
  
LamaConfig::LamaConfig()
{
    mCommunicationKind = lama::Matrix::SYNCHRONOUS;
    mComm              = lama::Communicator::get();
    mContextType       = memory::context::Host;
    mMaxIter           = lama::nIndex;
    mValueType         = common::scalar::DOUBLE;
    mLogLevel          = lama::LogLevel::convergenceHistory;
    mUseMetis          = false;
    mWeight            = 1.0f;
    mDevice            = 0;
}

LamaConfig::~LamaConfig()
{
}

bool LamaConfig::isNumber( const char* arg )
{
    int len = strlen( arg );

    for ( int i = 0; i < len; ++i )
    {
        if ( isdigit( arg[i] ) )
        {
            continue;
        }
        return false;
    }

    return true;
}

bool LamaConfig::isReal( const char* arg )
{
    int len = strlen( arg );

    for ( int i = 0; i < len; ++i )
    {
        if ( isdigit( arg[i] ) )
        {
            continue;
        }
        if ( arg[i] == '.' )
        {
            continue;
        }
        return false;
    }

    return true;
}

static void tokenize( std::vector<std::string>& tokens,
                      const std::string& str,
                      const std::string& delimiters = " ")
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of( delimiters, lastPos );

    while ( std::string::npos != pos || std::string::npos != lastPos )
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr( lastPos, pos - lastPos ) );
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of( delimiters, pos );
        // Find next "non-delimiter"
        pos = str.find_first_of( delimiters, lastPos );
    }
}

void LamaConfig::setArg( const char* arg )
{
    std::string val = arg;

    // check for multi option on a node, e.g. 'cpu,mic,gpu,gpu', each processor
    // on a node gets its own argument

    if ( strchr( arg, ',' ) )
    {
        std::vector<std::string> singleVals;

        tokenize( singleVals, val, "," );

        if ( singleVals.size() == mComm->getNodeSize() )
        {
            const std::string myVal = singleVals[ mComm->getNodeRank() ];

            setArg( myVal.c_str() );
        }
        else
        {
            std::cerr << val << " cannot be tokenized, #items = " << singleVals.size()
                      << " does not match node size = " << mComm->getNodeSize() << std::endl;
        }
        return;
    }

    // make upper string for more convenience, e.g. Host is same as host or HOST

    for ( std::string::iterator p = val.begin(); val.end() != p; ++p )
    {
            *p = toupper( *p );
    }

    if (   ( "CSR" == val ) || ( "ELL" == val ) || ( "COO" == val )
        || ( "DIA" == val ) || ( "JDS" == val ) )

    { 
        mMatrixFormat = val;
    }
    else if ( ( "HOST" == val ) || ( "CPU" == val ) )
    { 
        // Host does not require a device id

        mContextType = memory::context::Host;
    }
    else if ( ( "MIC" == val ) || ( "PHI" == val ) )
    { 
        mContextType = memory::context::MIC;
    }
    else if ( ( "CUDA" == val ) || ( "GPU" == val ) )
    { 
        mContextType = memory::context::CUDA;
    }
    else if ( "PINNED" == val )
    {
        // support fast memory transfer Host->CUDA

#ifdef USE_CUDA
        lama::CUDAHostContextManager::setAsCurrent( getContextPtr() );
#endif
    }
    else if ( "METIS" == val )
    {
        mUseMetis = true;
    }
    else if ( "SYNC" == val )
    {
        mCommunicationKind = lama::Matrix::SYNCHRONOUS;
    }
    else if ( "ASYNC" == val )
    {
        mCommunicationKind = lama::Matrix::ASYNCHRONOUS;
    }
    else if ( ( "FLOAT" == val ) || ( "SP" == val ) )
    {
        mValueType = common::scalar::FLOAT;
    }
    else if ( ( "DOUBLE" == val ) || ( "DP" == val ) )
    {
        mValueType = common::scalar::DOUBLE;
    }
    else if ( ( "LONGDOUBLE" == val ) || ( "LP" == val ) )
    {
        mValueType = common::scalar::LONG_DOUBLE;
    }
    else if ( ( "COMPLEX" == val ) || ( "CP" == val ) )
    {
        mValueType = common::scalar::COMPLEX;
    }
    else if ( ( "DOUBLECOMPLEX" == val ) || ( "COMPLEXDOUBLE" == val) || ( "ZP" == val ) )
    {
        mValueType = common::scalar::DOUBLE_COMPLEX;
    }
    else if ( ( "LONGCOMPLEX" == val ) || ( "COMPLEXLONG" == val) )
    {
        mValueType = common::scalar::LONG_DOUBLE_COMPLEX;
    }
    else if ( "TEXTURE" == val )
    {
        int replace = 1;
        setenv( "LAMA_CUDA_USE_TEXTURE", "1", replace );
    }
    else if ( "NOTEXTURE" == val )
    {
        int replace = 1;
        setenv( "LAMA_CUDA_USE_TEXTURE", "0", replace );
    }
    else if ( ( "SHAREDMEM" == val ) || ( "SM" == val ) )
    {
        int replace = 1;
        setenv( "LAMA_CUDA_USE_SHARED_MEM", "1", replace );
    }
    else if ( ( "NOSHAREDMEM" == val ) || ( "NOSM" == val ) )
    {
        int replace = 1;
        setenv( "LAMA_CUDA_USE_SHARED_MEM", "0", replace );
    }
    else if ( "LOG_HISTORY" == val ) 
    {
        mLogLevel = lama::LogLevel::convergenceHistory;
    }
    else if ( "LOG_SOLVER" == val ) 
    {
        mLogLevel = lama::LogLevel::solverInformation;
    }
    else if ( "LOG_AVANCED" == val ) 
    {
        mLogLevel = lama::LogLevel::advancedInformation;
    }
    else if ( "LOG_COMPLETE" == val ) 
    {
        mLogLevel = lama::LogLevel::completeInformation;
    }
    else if ( "LOG_NO" == val ) 
    {
        mLogLevel = lama::LogLevel::noLogging;
    }
    else if ( ( 'T' == val[0] ) && isNumber( val.c_str() + 1 ) )
    {
        int numThreads;
        int narg = sscanf( val.c_str() + 1, "%d", &numThreads );
        if ( narg > 0 )
        {
            omp_set_num_threads( numThreads );
        }
        else
        {
            std::cout << "Illegal for number of threads: " << arg << std::endl;
        }
    }
    else if ( ( 'D' == val[0] ) && isNumber( val.c_str() + 1 ) )
    {
        int deviceArg;

        int narg = sscanf( val.c_str() + 1, "%d", &deviceArg );
        if ( narg > 0 )
        {
            mDevice = deviceArg;
        }
        else
        {
            std::cout << "Illegal for number of device: " << arg << std::endl;
        }
    }
    else if ( ( 'B' == val[0] ) && isNumber( val.c_str() + 1 ) )
    {
        int numBlocks;
        int narg = sscanf( val.c_str() + 1, "%d", &numBlocks );
        if ( narg > 0 )
        {
            char envSetting[ 256 ];
            sprintf( envSetting, "LAMA_CUDA_BLOCK_SIZE=%d", numBlocks );
            int replace = 1;
            setenv( "LAMA_CUDA_BLOCK_SIZE", envSetting, replace );
            std::cout << "Environment setting: LAMA_CUDA_BLOCK_SIZE = " << envSetting << std::endl;
        }
        else
        {
            std::cout << "Illegal for block size on CUDA: " << arg << std::endl;
        }
    }
    else if ( ( 'W' == val[0] ) && isReal( val.c_str() + 1 ) )
    {
        float weight;

        int narg = sscanf( val.c_str() + 1, "%f", &weight );

        if ( narg > 0 && weight >= 0.0f )
        {
            mWeight = weight;
        }
        else
        {
            std::cout << "Illegal weight: " << arg << std::endl;
        }
    }
    else if ( isNumber( val.c_str() ) )
    {
        sscanf( val.c_str(), "%d", &mMaxIter );
    }
    else
    {
        std::cout << "Illegal argument: " << arg << std::endl;
    }
}

void LamaConfig::writeAt( std::ostream& stream ) const
{
    stream << "LAMA configuration" << std::endl;
    stream << "==================" << std::endl;
    stream << "Context           = " << getContext() << std::endl;
    stream << "Communicator      = " << *mComm << std::endl;
    stream << "Matrix format     = " << getFormat() << std::endl;
    stream << "CommKind          = " << mCommunicationKind << std::endl;
    stream << "ValueType         = " << mValueType << std::endl;
    stream << "#Threads/CPU      = " << omp_get_max_threads() << std::endl;
    stream << "weight            = " << mWeight << std::endl;

    if ( getenv( "LAMA_CUDA_USE_TEXTURE" ) )
    {
        stream << "useTexture(GPU)   = " << getenv( "LAMA_CUDA_USE_TEXTURE" ) << std::endl;
    }
    if ( getenv( "LAMA_CUDA_USE_SHARED_MEM" ) )
    {
        stream << "useSharedMem(GPU) = " << getenv( "LAMA_CUDA_USE_SHARED_MEM" ) << std::endl;
    }
    if ( getenv( "LAMA_CUDA_BLOCK_SIZE" ) )
    {
        stream << "BlockSize(GPU)    = " << getenv( "LAMA_CUDA_BLOCK_SIZE" ) << std::endl;
    }
    if ( hasMaxIter () )
    {
        stream << "Max iter          = " << mMaxIter << std::endl;
    }
}

template<typename ValueType>
lama::SparseMatrix<ValueType>* LamaConfig::createSparseMatrix()
{
    return createSparseMatrix<ValueType>( getFormat() );
}

const char* LamaConfig::getFormat( ) const
{
    if ( mMatrixFormat == "" )
    {
        // choose default format by context: Host -> CSR, CUDA -> ELL

        if ( mContextType == memory::context::CUDA )
        {
            return "ELL";
        }
        else
        {
            return "CSR";
        }
    }
    else
    {
        return mMatrixFormat.c_str();
    }
}

template<typename ValueType>
lama::SparseMatrix<ValueType>* LamaConfig::createSparseMatrix( const char* format )
{
    if ( strcmp( format, "CSR" ) == 0 )
    {
        return new lama::CSRSparseMatrix<ValueType>();
    }
    else if ( strcmp( format, "ELL" ) == 0 )
    {
        return new lama::ELLSparseMatrix<ValueType>();
    }
    else if ( strcmp( format, "JDS" ) == 0 )
    {
        return new lama::JDSSparseMatrix<ValueType>();
    }
    else if ( strcmp( format, "DIA" ) == 0 )
    {
        return new lama::DIASparseMatrix<ValueType>();
    }
    else if ( strcmp( format, "COO" ) == 0 )
    {
        return new lama::COOSparseMatrix<ValueType>();
    }
    else
    {
        std::string dformat = getFormat();

        std::cerr << "createSparseMatrix: " << format << " not supported, take " << dformat;

        return createSparseMatrix<ValueType>( dformat.c_str() );
    }
}

