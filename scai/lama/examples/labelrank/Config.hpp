/**
 * @file lama/examples/labelrank/Config.hpp
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
 * @brief Structure that contains configuration for label propagation
 * @author Thomas Brandes
 * @date 07.06.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/_Matrix.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/common/Printable.hpp>

#include <cstring>
#include <locale>

/** Parameters to define LAMA behavior */

class Config : public scai::common::Printable
{

public:

    Config()
    {
        // overlap communication with local computation
        mCommunicationKind = scai::lama::SyncKind::SYNCHRONOUS;
        mComm              = scai::dmemo::Communicator::getCommunicatorPtr();
        mContext           = scai::hmemo::Context::getHostPtr();
        mMaxIters          = 1000;
    }

    ~Config()
    {
        // give up ownership for communicator and context
        mComm.reset();
        mContext.reset();
    }

    void setArg( const char* arg )
    {
        std::string val = arg;
        std::locale loc;

        // make upper string for more convenience, e.g. Host is same as host or HOST

        for ( std::string::iterator p = val.begin(); val.end() != p; ++p )
        {
            *p = std::toupper( *p, loc );
        }

        // make it upper case

        if (   ( "CSR" == val ) || ( "Dense" == val ) )
        {
            mMatrixFormat = val;
        }
        else if ( "HOST" == val )
        {
            mContext = scai::hmemo::Context::getContextPtr( scai::common::ContextType::Host );
        }
        else if ( ( "CUDA" == val ) || ( "GPU" == val ) )
        {
            // int device = mComm->getNodeRank();
            int device = 0;
            mContext = scai::hmemo::Context::getContextPtr( scai::common::ContextType::CUDA, device );
        }
        else if ( "SYNC" == val )
        {
            mCommunicationKind = scai::lama::SyncKind::SYNCHRONOUS;
        }
        else if ( "ASYNC" == val )
        {
            mCommunicationKind = scai::lama::SyncKind::ASYNC_LOCAL;
        }
        else if ( isNumber( val.c_str() ) )
        {
            sscanf( val.c_str(), "%d", &mMaxIters );
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
            if ( mContext->getType() == scai::common::ContextType::CUDA )
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

    scai::hmemo::ContextPtr getContextPtr() const
    {
        return mContext;
    }

    const scai::hmemo::Context& getContext() const
    {
        return *mContext;
    }

    scai::dmemo::CommunicatorPtr getCommunicatorPtr() const
    {
        return mComm;
    }

    const scai::dmemo::Communicator& getCommunicator() const
    {
        return *mComm;
    }

    std::string             mMatrixFormat;
    scai::hmemo::ContextPtr mContext;
    scai::lama::SyncKind    mCommunicationKind;
    int                     mMaxIters;

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

    scai::dmemo::CommunicatorPtr  mComm;

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

