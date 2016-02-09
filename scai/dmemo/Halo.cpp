/**
 * @file Halo.cpp
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
 * @brief Halo.cpp
 * @author Thomas Brandes
 * @date 23.02.2011
 * @since 1.0.0
 */

// hpp
#include <scai/dmemo/Halo.hpp>

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( Halo::logger, "Halo" )

Halo::Halo()
{
}

Halo::Halo( const Halo& other )
{
    operator=( other );
}

Halo::~Halo()
{
}

void Halo::clear()
{
    mRequiredPlan.clear();
    mProvidesPlan.clear();

    mRequiredIndexes.clear();
    mProvidesIndexes.clear();

    mGlobal2Halo.clear();
}

void Halo::purge()
{
    mRequiredPlan.purge();
    mProvidesPlan.purge();

    mRequiredIndexes.purge();
    mProvidesIndexes.purge();

    // free memory of map by reallocation

    std::map<IndexType,IndexType>().swap( mGlobal2Halo );
}

Halo& Halo::operator=( const Halo& other )
{
    if( this != &other )
    {
        SCAI_LOG_DEBUG( logger, "make deep copy of Halo" )

        mRequiredPlan = other.mRequiredPlan;
        mProvidesPlan = other.mProvidesPlan;

        mRequiredIndexes = other.mRequiredIndexes;
        mProvidesIndexes = other.mProvidesIndexes;

        mGlobal2Halo = other.mGlobal2Halo;
    }

    return *this;
}

/* ---------------------------------------------------------------------- */

void Halo::writeAt( std::ostream& stream ) const
{
    // write info this object

    stream << "Halo( size = " << getHaloSize() << ", required plan = " << mRequiredPlan << ", provides plan = "
           << mProvidesPlan << ")";
}

} /* end namespace dmemo */

} /* end namespace scai */
