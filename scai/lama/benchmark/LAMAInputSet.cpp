/**
 * @file LAMAInputSet.cpp
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
 * @brief LAMAInputSet.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 12.09.2011
 */

#include <scai/lama/benchmark/LAMAInputSet.hpp>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( LAMAInputSet::logger, "InputSet.LAMAInputSet" );

LAMAInputSet::LAMAInputSet() : InputSet( "", "LAMAInputSet" )
{
    mAlpha = 1.0;
    mBeta  = 1.0;

    SCAI_LOG_INFO( logger, "LAMA input set created, all elems are NULL" );
}

LAMAInputSet::~LAMAInputSet()
{
    SCAI_LOG_INFO( logger, "~LAMAInputSet, input data is freed via unique pointers" );
}

void LAMAInputSet::writeAt( std::ostream& stream ) const
{
    stream << "LAMAInputSet( alpha = " << mAlpha << ", beta = " << mBeta;

    if ( mA.get() )
    {
        stream << ", A = " << *mA;
    }
    else
    {
        stream << ", A = -"; 
    }

    if ( mX.get() )
    {
        stream << ", X = " << *mX;
    }
    else
    {
        stream << ", X = -"; 
    }

    if ( mY.get() )
    {
        stream << ", Y = " << *mY;
    }
    else
    {
        stream << ", Y = -";
    }

    stream << " )";
}

double LAMAInputSet::getAlpha() const
{
    return mAlpha;
}

double LAMAInputSet::getBeta() const
{
    return mBeta;
}

const lama::DenseVector<double>& LAMAInputSet::getX() const
{
    SCAI_ASSERT_ERROR( mX.get(), "mX not set before" )

    return *mX;
}

const lama::DenseVector<double>& LAMAInputSet::getY() const
{
    SCAI_ASSERT_ERROR( mY.get(), "mY not set before" )

    return *mY;
}

const lama::CSRSparseMatrix<double>& LAMAInputSet::getA() const
{
    SCAI_ASSERT_ERROR( mA.get(), "mA not set before" )

    return *mA;
}

}

}
