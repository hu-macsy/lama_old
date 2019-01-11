/**
 * @file Create.cpp
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
 * @brief Demo program for the Factory of HArray.
 * @author Thomas Brandes, Lauretta Schubert
 * @date 18.04.2012
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/throw.hpp>

#include <iostream>
#include <memory>

using namespace scai;
using namespace hmemo;

SCAI_LOG_DEF_LOGGER( logger, "CreateTest" )

// Template instantiation of LAMArray

template class scai::hmemo::HArray<double>;

int main()
{
    SCAI_LOG_THREAD( "Main" )
    ContextPtr contextPtr = Context::getHostPtr();
    static IndexType N =  100;
    HArray<float> lamaArray ( N, 1.0 );
    std::shared_ptr<HArray<float> > lamaArray1( lamaArray.newArray() );
    *lamaArray1 = lamaArray;
    ReadAccess<float> read( lamaArray, contextPtr );
    ReadAccess<float> read1( *lamaArray1, contextPtr );
    const float* data = read.get();
    const float* data1 = read1.get();

    for ( IndexType i = 0; i < N; ++i )
    {
        SCAI_ASSERT_EQUAL( data[i], data1[i], "" )
    }

    std::cout << "Create finished" << std::endl;
    std::shared_ptr<_HArray> lamaArray2( _HArray::create( scai::common::ScalarType::FLOAT ) );
    std::cout << "lamaArray2 = " << *lamaArray2 << std::endl;
    std::shared_ptr<_HArray> lamaArray3( _HArray::create( scai::common::ScalarType::DOUBLE ) );
    std::cout << "lamaArray3 = " << *lamaArray3 << std::endl;
}
