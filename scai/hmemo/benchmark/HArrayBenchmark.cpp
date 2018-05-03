/**
 * @file HArrayBenchmark.cpp
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
 * @brief Implementation of methods for benchmark with HArray
 * @author Thomas Brandes
 * @date 21.09.2017
 */

#include <scai/hmemo/benchmark/HArrayBenchmark.hpp>

#include <scai/benchmark/Parser.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>

using std::string;

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( HArrayBenchmark::logger, "Benchmark.HArrayBenchmark" );

HArrayBenchmark::HArrayBenchmark( const string& argument ) :

    Benchmark( "", getCreateId() ),
    mArgument( argument )
{
    if ( argument == "" )
    {
        mArgument = "Host, Host";
    }

    std::vector<string> argTokens;

    tokenize( argTokens, mArgument, " ,:" );

    SCAI_ASSERT_EQ_ERROR( argTokens.size(), 2, "HArrayBenchmark( context1, context2 ), two args expected" )

    common::ContextType ctx1 = common::str2ContextType( argTokens[0].c_str() );
    common::ContextType ctx2 = common::str2ContextType( argTokens[1].c_str() );

    if ( ctx1 == common::ContextType::MaxContext )
    {
        COMMON_THROWEXCEPTION( argTokens[0] << " does not specify a context" )
    }

    if ( ctx2 == common::ContextType::MaxContext )
    {
        COMMON_THROWEXCEPTION( argTokens[1] << " does not specify a context" )
    }

    mTargetContext = Context::getContextPtr( ctx1 );
    mSourceContext = Context::getContextPtr( ctx2 );

    mArgument = std::string( contextType2str( ctx1 ) ) + ", " + std::string( contextType2str( ctx2 ) );

    mName = getCreateId() + "( " +  mArgument + " )";

    SCAI_LOG_INFO( logger, "Benchmark " << mName << " created" )
}

HArrayBenchmark::~HArrayBenchmark()
{
}

void HArrayBenchmark::writeAt( std::ostream& stream ) const
{
    stream << "HArrayBenchmark( target context = " << *mTargetContext
           << ", source context = " << *mSourceContext << " )";
}

void HArrayBenchmark::initialize()
{
    mInputSet.reset( benchmark::InputSet::createWithArgument( mInputSetId ) );
  
    mHArrayInputSet = dynamic_cast<HArrayInputSet*>( mInputSet.get() );

    SCAI_ASSERT_ERROR( mHArrayInputSet, "Not a HArrayInputSet: " << *mInputSet )

    mArray = &mHArrayInputSet->getArray();

    SCAI_LOG_INFO( logger, "Initialized: benchmark = " << *this << ", input set = " << *mHArrayInputSet )
}

void HArrayBenchmark::setUp()
{
    // Allocate the memory
}

void HArrayBenchmark::execute()
{
    IndexType N = mArray->size();

    SCAI_LOG_TRACE( logger, "execute, size = " << N )

    // write only access on source context invalidates data at target context
    {   
        WriteOnlyAccess<double> write( *mArray, mSourceContext, N );
    }   
    // read access transfers data from context1 to context2
    {   
        ReadAccess<double> read( *mArray, mTargetContext );
    }
}

void HArrayBenchmark::tearDown()
{
}

void HArrayBenchmark::shutdown()
{
}

CounterType HArrayBenchmark::getNumFloatingPointOperations() const
{
    return 0;
}

CounterType HArrayBenchmark::getProcessedBytes() const
{
    return sizeof( double ) * mArray->size();
}

}

}
