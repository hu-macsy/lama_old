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

#include <lama.hpp>

#include <lama/Context.hpp>
#include <lama/Printable.hpp>
#include <lama/matrix/Matrix.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <cstring>

/** Class that handles commonly used configuration values for LAMA
 *
 *  Very important: DO NOT USE IT for a global variable
 *
 *  \code
 *      LAMAConfig lamaconf;    // constructur might fail 
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

    LamaConfig()
    {
        mCommunicationKind = lama::Matrix::SYNCHRONOUS;
        mComm              = lama::CommunicatorFactory::get();
        mContext           = lama::ContextFactory::getContext( lama::Context::Host );
    }

    ~LamaConfig()
    {
        // due to shared pointer destructor gives up ownership of context, communicator
    }

    void setArg( const char* arg )
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
            // int device = mComm->getNodeRank();
            int device = 0;
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
        else
        {
            std::cout << "Illegal argument: " << arg << std::endl;
        }
    }

    const char* getFormat( ) const
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

    lama::ContextPtr getContextPtr() const
    {
        return mContext;
    }

    lama::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const lama::Communicator& getCommunicator() const
    {
        return *mComm;
    }

    std::string            mMatrixFormat;
    lama::ContextPtr       mContext;
    lama::Matrix::SyncKind mCommunicationKind;

    void writeAt( std::ostream& stream ) const
    {
        stream << "LAMA configuration" << std::endl;
        stream << "==================" << std::endl;
        stream << "Context       = " << *mContext << std::endl;
        stream << "Communicator  = " << *mComm << std::endl;
        stream << "Matrix format = " << getFormat() << std::endl;
        stream << "CommKind      = " << mCommunicationKind << std::endl;
    }

private:

    lama::CommunicatorPtr  mComm;
};

