/**
 * @file HArrayInputSet.cpp
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
 
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of methods of class HArrayInputSet.
 * @author Thomas Brandes
 * @date 18.09.2017
 */

#include <scai/hmemo/benchmark/HArrayInputSet.hpp>

#include <scai/benchmark/Parser.hpp>

#include <string>

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( HArrayInputSet::logger, "InputSet.HArrayInputSet" );

HArrayInputSet::HArrayInputSet( const std::string& argument ) : 

    InputSet( "", "HArrayInputSet" ),
    mArgument( argument )
{
    if ( argument == "" )
    {
        mArgument = "Host 1000000";
    }

    std::vector<std::string> argTokens;

    tokenize( argTokens, mArgument, " ,:" );

    SCAI_ASSERT_EQ_ERROR( argTokens.size(), 2, "HArrayInputSet( memory, size ), two args expected" )

    common::ContextType ctx = common::str2ContextType( argTokens[0].c_str() );

    if ( ctx == common::ContextType::MaxContext )
    {
        COMMON_THROWEXCEPTION( argTokens[0] << " does not specify a context" )
    }

    mContext = Context::getContextPtr( ctx );

    mSize    = static_cast<IndexType>( strtol( argTokens[1].c_str(), NULL, 10 ) );

    std::ostringstream args; 
  
    args << ctx << ", " << mSize;

    mArgument = args.str();

    resetName( getCreateId() + "( " +  mArgument + " )" );

    mArray.reset( new HArray<double>( mSize, 0.0, mContext ) );

    SCAI_LOG_INFO( logger, "InputSet " << getName() << " created, array = " << *mArray )
}

HArrayInputSet::~HArrayInputSet()
{
    SCAI_LOG_INFO( logger, "~HArrayInputSet" );
}

void HArrayInputSet::writeAt( std::ostream& stream ) const
{
    // *mArray is available over the whole lifetime of the input set

    stream << "HArrayInputSet( array = " << *mArray << " )";
}

HArray<double>& HArrayInputSet::getArray()
{
    SCAI_ASSERT_ERROR( mArray.get(), "no array allocated" );
    return *mArray;
}

}

}
