/**
 * @file tracing/examples/TraceTest.cpp
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
 * @brief ToDo: Missing description in ./tracing/examples/TraceTest.cpp
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#include <scai/tracing.hpp>
#include <scai/logging.hpp>
#include <scai/common/Walltime.hpp>

void subA()
{
    SCAI_REGION( "A" )
    scai::common::Walltime::sleep( 1000 );
}

void subB()
{
    SCAI_REGION( "B" )
    subA();
    scai::common::Walltime::sleep( 2000 );
}

int main()
{
    SCAI_LOG_THREAD( "master" )
    SCAI_REGION( "main" )
    subA();
    subB();
    scai::common::Walltime::sleep( 3000 );
}

