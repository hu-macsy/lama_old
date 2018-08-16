/**
 * @file macros/tracing.hpp
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
 * @brief Definition of the tracing macros.
 * @author Lauretta Schubert, Thomas Brandes
 * @date 12.08.2011
 */

#pragma once

#if defined( SCAI_TRACE_ON )

#include <scai/tracing/TraceRegionRecord.hpp>
#include <scai/tracing/TraceConfig.hpp>

#define SCAI_REGION( name ) scai::tracing::ScopedTraceRecord SCAI_Trc__( name, __FILE__, __LINE__ );
#define SCAI_REGION_N( name, n ) scai::tracing::ScopedTraceRecord SCAI_Trc__( name, n, __FILE__, __LINE__ );
#define SCAI_REGION_START( name ) scai::tracing::TraceRegionRecord::start( name, __FILE__, __LINE__ );
#define SCAI_REGION_END( name ) scai::tracing::TraceRegionRecord::stop( name );
#define SCAI_TRACE_SCOPE( flag ) scai::tracing::TraceConfig::TraceScope SCAI_Scp__( flag );

#elif defined( SCAI_TRACE_OFF )

#define SCAI_REGION( name )
#define SCAI_REGION_START( name )
#define SCAI_REGION_END( name )
#define SCAI_REGION_N( name, n )
#define SCAI_TRACE_SCOPE( flag )

#else

//Macro SCAI_REGION should also be defined in case of error for convenience with Eclipse.
#define SCAI_REGION( name )
#define SCAI_REGION_START( name )
#define SCAI_REGION_END( name )
#define SCAI_REGION_N( name, n )
#define SCAI_TRACE_SCOPE( flag )
// turned off for master branch
// #pragma message("Must define SCAI_TRACE_xxx with xxx = ON, or OFF")

#endif
