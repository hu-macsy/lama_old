/**
 * @file GPIMemManager.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Management of GASPI segment data
 * @author Thomas Brandes
 * @date 09.05.2014
 * @since 1.1.0
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

