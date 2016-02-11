/**
 * @file Distributed.cpp
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
 * @brief Distributed.cpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */

// hpp
#include <scai/dmemo/Distributed.hpp>

namespace scai
{

namespace dmemo
{

Distributed::Distributed( DistributionPtr distribution )

    : mDistribution( distribution )
{
    if ( !distribution )
    {
        COMMON_THROWEXCEPTION( "Distributed object must not have NULL distribution" )
    }
}

Distributed::Distributed( const Distributed& other )

    : mDistribution( other.mDistribution )
{
    // copy shared pointer is okay, mDistribution can never be NULL
}

Distributed::~Distributed()
{
}

void Distributed::swap( Distributed& other )
{
    mDistribution.swap( other.mDistribution );
}

void Distributed::setDistributionPtr( DistributionPtr distributionPtr )
{
    mDistribution = distributionPtr;
}

void Distributed::buildCSRGraph( IndexType ia[], IndexType[], IndexType vwgt[], const IndexType* ) const
{
    IndexType size = mDistribution->getLocalSize();

    for ( IndexType i = 0; i < size; ++i )
    {
        ia[i] = 0;
        vwgt[i] = 0;
    }
}

IndexType Distributed::getCSRGraphSize() const
{
    return 0;
}

} /* end namespace dmemo */

} /* end namespace scai */
