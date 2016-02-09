/**
 * @file AMGSetup.cpp
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
 * @brief AMGSetup.cpp
 * @author Jiri Kraus
 * @date 28.10.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/AMGSetup.hpp>

namespace scai
{

namespace solver
{

AMGSetup::AMGSetup()
    : mHostOnlyLevel( std::numeric_limits<IndexType>::max() ), mHostOnlyVars( 0 ), mReplicatedLevel(
          std::numeric_limits<IndexType>::max() )
{
}

AMGSetup::~AMGSetup()
{
}

void AMGSetup::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    mHostOnlyLevel = hostOnlyLevel;
}

void AMGSetup::setHostOnlyVars( IndexType hostOnlyVars )
{
    mHostOnlyLevel = hostOnlyVars;
}

void AMGSetup::setReplicatedLevel( IndexType replicatedLevel )
{
    mReplicatedLevel = replicatedLevel;
}

void AMGSetup::writeAt( std::ostream& stream ) const
{
    stream << "AMGSetup( ... )";
}

} /* end namespace solver */

} /* end namespace scai */
