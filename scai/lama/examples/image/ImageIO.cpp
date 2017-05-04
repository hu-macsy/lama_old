/**
 * @file ImageIO.cpp
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
 * @brief Implementation of routines to read/write image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/GridVector.hpp>

#include <scai/lama/examples/image/ImageIO.hpp>

namespace scai
{

using namespace utilskernel;

namespace lama
{

SCAI_LOG_DEF_LOGGER( ImageIO::logger, "ImageIO" )

template<typename ValueType>
void ImageIO::read( GridVector<ValueType>& imageData, const std::string& inputFileName )
{
    IndexType width = 10;
    IndexType height = 10;

    common::Grid3D grid( width, height, 3 );

    LArray<ValueType> x( grid.size() );

    imageData.swap( x, grid );

    SCAI_LOG_ERROR( logger, "cannot read data from file " << inputFileName )
}

template<typename ValueType>
void ImageIO::write( const GridVector<ValueType>& imageData, const std::string& outputFileName )
{
    SCAI_LOG_ERROR( logger, "cannot write data to file " << outputFileName )

    const LArray<ValueType>& x = imageData.getLocalValues();

    COMMON_THROWEXCEPTION( "write image not yet supported, data = " << x )
}

// instantiate methods

template
void ImageIO::write( const GridVector<float>&, const std::string& );

template
void ImageIO::read( GridVector<float>&, const std::string& );

} /* end namespace lama */

} /* end namespace scai */
