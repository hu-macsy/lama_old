/**
 * @file GPIMemManager.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Management of GASPI segment data
 * @author Thomas Brandes
 * @date 09.05.2014
 */

#pragma once

// for GASPI data types

#include <GASPI.h>

// for dll_import
#include <scai/common/config.hpp>

#include <vector>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace dmemo
{

/** Class with static methods for memory management */

class COMMON_DLL_IMPORTEXPORT GPIMemManager
{

public:

    /** Get unused segment data for user application. 
     *
     *  @param[out] id is the identification of GASPI segment
     *  @param[out] ptr is pointer to the segment data
     *  @param[out] offset offset within the used segment
     *  @param[in]  size is the number of bytes required for the segment data
     *
     *  This routine might assign data of the same segment for two following calls so offset is 
     *  used to distinguish between the different data.
     *
     *  NOTE: currently offset is always 0 as otherwise consistency of the argument size between
     *        all GASPI processes must be guaranteed.
     */

    static void getSegmentData( gaspi_segment_id_t& id, gaspi_pointer_t& ptr, gaspi_offset_t& offset, const int size );

    /** This method checks for a given pointer in which segment it is available.
     *  
     *  @param[out] id is the identification of GASPI segment
     *  @param[out] offset offset within the used segment
     *  @param[in]  ptr is pointer to be searched for
     *
     *  @return true if ptr belongs to a allocated and used GASPI segment
     *
     *  If this method returns false, id and offset are undefined
     */
    static bool findSegment( gaspi_segment_id_t& id, gaspi_offset_t& offset, const gaspi_pointer_t ptr );

    /** Release used segment data from application. */

    static void releaseSegmentData( const gaspi_segment_id_t id, const gaspi_offset_t offset );

    /** This routine frees all segment data. It should not throw an exception. */

    static void freeAll();

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace dmemo

} // namespace scai

