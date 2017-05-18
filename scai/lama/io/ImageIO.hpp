/**
 * @file ImageIO.hpp
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

    /** Read in a grid vector from an input file.
     *
     *  @param[in] gridData will contain the grid vector data.
     *  @param[in] inputFileName name of the input file.
     *
     */
    template<typename ValueType>
    static void read( GridVector<ValueType>& gridData, const std::string& inputFileName );

    /** Write the pixel data to an output file.
     *
     *  Supported formats: *.png
     */
    template<typename ValueType>
    static void write( const GridVector<ValueType>& gridData, const std::string& outputFileName );

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

/** Metaprogramming structure to call a routine for each type in a typelist 
 *
 *  This struct should be used by all derived classes.
 */

template<class Derived, typename TList>
struct ImageIOWrapper;

/*
 * Termination
 */
template<class Derived>
struct ImageIOWrapper<Derived, common::mepr::NullType>
{
    static void write( Derived&, const hmemo::_HArray& data, const common::Grid&, const std::string& )
    {
        COMMON_THROWEXCEPTION( "write " << data << " unsupported, unknown type." )
    }

    static void read( Derived&, hmemo::_HArray& data, common::Grid&, const std::string& )
    {
        COMMON_THROWEXCEPTION( "read " << data << " unsupported, unknown type." )
    }
};

template<class Derived, typename ValueType, typename TailTypes>
struct ImageIOWrapper<Derived, common::mepr::TypeList<ValueType, TailTypes> >
{
    static void write( Derived& io, const hmemo::_HArray& data, const common::Grid& grid, const std::string& fileName )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeImpl( reinterpret_cast<const hmemo::HArray<ValueType>& >( data ), grid, fileName );
        }
        else
        {
            ImageIOWrapper<Derived, TailTypes>::write( io, data, grid, fileName );
        }
    }

    static void read( Derived& io, hmemo::_HArray& data, common::Grid& grid, const std::string& fileName )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readImpl( reinterpret_cast<hmemo::HArray<ValueType>& >( data ), grid, fileName );
        }
        else
        {
            ImageIOWrapper<Derived, TailTypes>::read( io, data, grid, fileName );
        }
    }
};


} /* end namespace lama */

} /* end namespace scai */
