/**
 * @file NoDistribution.cpp
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
 * @brief Implementation of methods for class NoDistribution.
 * @author Thomas Brandes
 * @date 14.03.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/distribution/NoDistribution.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>

// std
#include <fstream>

#define MASTER 0

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( NoDistribution::logger, "Distribution.NoDistribution" )

NoDistribution::NoDistribution( const IndexType globalSize )
    : Distribution( globalSize )
{
}

NoDistribution::~NoDistribution()
{
}

bool NoDistribution::isLocal( const IndexType /* index */) const
{
    return true;
}

IndexType NoDistribution::getLocalSize() const
{
    return mGlobalSize;
}

IndexType NoDistribution::local2global( const IndexType localIndex ) const
{
    return localIndex;
}

IndexType NoDistribution::global2local( const IndexType globalIndex ) const
{
    return globalIndex;
}

bool NoDistribution::isEqual( const Distribution& other ) const
{
    return typeid( *this ) == typeid( other ) && getGlobalSize() == other.getGlobalSize();
}

void NoDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "NoDistribution( size = " << mGlobalSize << " )";
}

void NoDistribution::printDistributionVector( std::string name ) const
{
    if( mCommunicator->getRank() == MASTER ) // process 0 ist MASTER process
    {
        std::ofstream file;
        file.open( ( name + ".part" ).c_str() );
        // print row - partition mapping
        file << "No Distribution: all rows are available on all processes." << std::endl;
        file.close();
    }
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string NoDistribution::createValue()
{
    return "NO";
}

Distribution* NoDistribution::create( const DistributionArguments arg )
{
    // Note: weight argument is not used here
    //       same is true for matrix, commonunicationPtr

    return new NoDistribution( arg.globalSize );
}

} /* end namespace lama */

} /* end namespace scai */
