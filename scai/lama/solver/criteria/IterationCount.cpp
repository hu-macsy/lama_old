/**
 * @file IterationCount.cpp
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
 * @brief IterationCount.cpp
 * @author Kai Buschulte
 * @date 21.07.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/criteria/IterationCount.hpp>

// others
#include <scai/lama/solver/IterativeSolver.hpp>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( IterationCount::logger, "Criterion.IterationCount" );

IterationCount::IterationCount()
    : Criterion(), mIterationExtrema( 1 )
{
}

IterationCount::IterationCount( const IndexType iterationExtrema )
    : Criterion(), mIterationExtrema( iterationExtrema )
{
    SCAI_LOG_DEBUG( logger, "Creating IterationCount with " << iterationExtrema );
}

IterationCount::IterationCount( const IterationCount &other )
    : Criterion( other ), mIterationExtrema( other.mIterationExtrema )

{
}

IterationCount::~IterationCount()
{
}

bool IterationCount::isSatisfied( const scai::lama::IterativeSolver& solver )
{
    SCAI_LOG_INFO( logger,
                   "Iteration Extrema = " << mIterationExtrema << ", Iteration Count = " << solver.getIterationCount() );

    return solver.getIterationCount() >= mIterationExtrema;
}

IndexType IterationCount::getIterationExtrema() const
{
    return mIterationExtrema;
}

void IterationCount::setIterationExtrema( IndexType iterationExtrema )
{
    mIterationExtrema = iterationExtrema;
}

void IterationCount::writeAt( std::ostream& stream ) const
{
    stream << "ItCount<" << getIterationExtrema() << ">";
}

} /* end namespace lama */

} /* end namespace scai */
