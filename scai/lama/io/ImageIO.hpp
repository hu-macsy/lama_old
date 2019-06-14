/**
 * @file lama/io/ImageIO.hpp
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
 * @brief Definition of routines to read/write image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/logging.hpp>
#include <scai/lama/GridVector.hpp>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT ImageIO
{

public:

    /** Write a two-dimensional array as scaled image. */

    template<typename ValueType>
    static void writeSC( const GridVector<ValueType>& arrayData, const std::string& outputFileName );

    /** Write a two-dimensional array as scaled image with given minimal and maximal value
     *
     *  @param[in] arrayData must be a two-dimensional grid
     *  @param[in] minVal is used as minimal value for color scaling
     *  @param[in] maxVal is used as maximal value for color scaling
     *  @param[in] outputFileName is the name of output file, suffix decides about format
     */
    template<typename ValueType>
    static void writeSC( const GridVector<ValueType>& arrayData, const ValueType minVal, const ValueType maxVal, const std::string& outputFileName );

    /** Implementation of virtual routine ImageIO::write for this format */

    virtual void read( hmemo::_HArray& data, common::Grid& grid, const std::string& outputFileName ) = 0;

    /** Implementation of virtual routine ImageIO::write for this format */

    virtual void write( const hmemo::_HArray& data, const common::Grid& grid, const std::string& outputFileName ) = 0;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
