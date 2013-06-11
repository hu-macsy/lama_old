/**
 * @file LAMAConfig.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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

#ifndef LAMA_CONFIG_HPP_
#define LAMA_CONFIG_HPP_

#include <lama.hpp>

#include <lama/Context.hpp>
#include <lama/Printable.hpp>
#include <lama/matrix/all.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <cstring>

/** Class that handles commonly used configuration values for LAMA
 *
 *  Very important: DO NOT USE THIS CLASS for a GLOBAL variable
 *
 *  ( so its constructor might be called before LAMA factories are available )
 *
 *  \code
 *      LAMAConfig lamaconf;    // constructur might fail 
 *
 *      int main () 
 *      {
 *          LAMAConfig lamaconf();   // that is the right use
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
    lama::SparseMatrix<ValueType>* createSparseMatrix()
    {
        return createSparseMatrix<ValueType>( getFormat() );
    }

    template<typename ValueType>
    lama::SparseMatrix<ValueType>* createSparseMatrix( const char* format );

    lama::ContextPtr getContextPtr() const
    {
        return mContext;
    }

    const lama::Context& getContext() const
    {
        return *mContext;
    }

    lama::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const lama::Communicator& getCommunicator() const
    {
        return *mComm;
    }

    std::string              mMatrixFormat;
    lama::ContextPtr         mContext;
    lama::Matrix::SyncKind   mCommunicationKind;
    lama::Scalar::ScalarType mValueType;          // value type to use

    void writeAt( std::ostream& stream ) const
    {
        stream << "LAMA configuration" << std::endl;
        stream << "==================" << std::endl;
        stream << "Context       = " << *mContext << std::endl;
        stream << "Communicator  = " << *mComm << std::endl;
        stream << "Matrix format = " << getFormat() << std::endl;
        stream << "CommKind      = " << mCommunicationKind << std::endl;
        stream << "ValueType     = " << mValueType << std::endl;
        if ( hasMaxIter () )
        {
            stream << "Max iter      = " << mMaxIter << std::endl;
        }
    }

    bool hasMaxIter() const
    { 
        return mMaxIter != lama::nIndex; 
    }

    lama::IndexType getMaxIter() const
    { 
        return mMaxIter; 
    }

private:

    lama::CommunicatorPtr  mComm;

    lama::IndexType        mMaxIter;

    inline bool isNumber( const char* arg )
    {
        int len = strlen( arg );

        for ( int i = 0; i < len; ++i )
        {
            if ( !isdigit( arg[i] ) )
            {
                return false;
            }
        }

        return true;
    }
};

/* ---------------------------------------------------------------------------- */
  
LamaConfig::LamaConfig()
{
    mCommunicationKind = lama::Matrix::SYNCHRONOUS;
    mComm              = lama::CommunicatorFactory::get();
    mContext           = lama::ContextFactory::getContext( lama::Context::Host );
    mMaxIter           = lama::nIndex;
    mValueType         = lama::Scalar::DOUBLE;
}

LamaConfig::~LamaConfig()
{
}

void LamaConfig::setArg( const char* arg )
{
    std::string val = arg;

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
    else if ( "HOST" == val )
    { 
        mContext = lama::ContextFactory::getContext( lama::Context::Host );
    }
    else if ( ( "CUDA" == val ) || ( "GPU" == val ) )
    { 
        // int device = 0;
        // adapt it for your node configuration

        int device = mComm->getNodeRank();

        mContext = lama::ContextFactory::getContext( lama::Context::CUDA, device );
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
        mValueType = lama::Scalar::FLOAT;
    }
    else if ( ( "DOUBLE" == val ) || ( "DP" == val ) )
    {
        mValueType = lama::Scalar::DOUBLE;
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

const char* LamaConfig::getFormat( ) const
{
    if ( mMatrixFormat == "" )
    {
        // choose default format by context: Host -> CSR, CUDA -> ELL

        if ( mContext->getType() == lama::Context::CUDA )
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


#endif
