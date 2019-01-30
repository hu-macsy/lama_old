/**
 * @file IOWrapper.hpp
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
 * @brief IOwrapper struct as metaprogramming to call typed IO routines
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#pragma once

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

#include <cstdio>

namespace scai
{

namespace lama
{

/* --------------------------------------------------------------------------------- */
/*  Metaprogramming for different Value Types                                        */
/* --------------------------------------------------------------------------------- */

/** Metaprogramming structure to call a routine for each type in a typelist 
 *
 *  \code
 *      IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( *this, array );
 *  \endcode
 */
template<class IOClass, typename TList>
struct IOWrapper;

/*
 * Termination
 */
template<class IOClass>
struct IOWrapper<IOClass, common::mepr::NullType>
{
    static void writeStorageImpl( IOClass&, const _MatrixStorage& storage )
    {
        COMMON_THROWEXCEPTION( "writeStorage " << storage << " unsupported, unknown type." )
    }

    static void readStorageImpl( IOClass&, _MatrixStorage& storage )
    {
        COMMON_THROWEXCEPTION( "readStorage " << storage << " unsupported, unknown type." )
    }

    static void writeArrayImpl( IOClass&, const hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void writeSparseImpl( IOClass&, const IndexType, const hmemo::HArray<IndexType>&, const hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void readArrayImpl( IOClass&, hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "readArray " << array << " unsupported, unknown type." )
    }

    static void readSparseImpl( IOClass&, IndexType&, hmemo::HArray<IndexType>&, hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "readSparse unsupported, " << array.getValueType() << " not in SCAI_ARRAY_TYPES"  )
    }

    static void writeGridImpl( IOClass&, const hmemo::_HArray& data, const common::Grid& )
    {
        COMMON_THROWEXCEPTION( "write " << data << " unsupported, unknown type." )
    }

    static void readGridImpl( IOClass&, hmemo::_HArray& data, common::Grid& )
    {
        COMMON_THROWEXCEPTION( "read " << data << " unsupported, unknown type." )
    }
};

/*
 * Step n
 */
template<class IOClass, typename ValueType, typename TailTypes>
struct IOWrapper<IOClass, common::mepr::TypeList<ValueType, TailTypes> >
{
    static void writeStorageImpl( IOClass& io, const _MatrixStorage& storage )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeStorageImpl( static_cast<const MatrixStorage<ValueType>& >( storage ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeStorageImpl( io, storage );
        }
    }

    static void readStorageImpl( IOClass& io, _MatrixStorage& storage )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readStorageImpl( static_cast< MatrixStorage<ValueType>& >( storage ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readStorageImpl( io, storage );
        }
    }

    static void writeArrayImpl( IOClass& io, const hmemo::_HArray& array )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeArrayImpl( static_cast<const hmemo::HArray<ValueType>& >( array ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeArrayImpl( io, array );
        }
    }

    static void writeSparseImpl( 
        IOClass& io, 
        const IndexType size, 
        const hmemo::HArray<IndexType>& indexes, 
        const hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeSparseImpl( size, indexes, static_cast<const hmemo::HArray<ValueType>& >( values ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeSparseImpl( io, size, indexes, values );
        }
    }

    static void readArrayImpl( IOClass& io, hmemo::_HArray& array )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readArrayImpl( static_cast< hmemo::HArray<ValueType>& >( array ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readArrayImpl( io, array );
        }
    }

    static void readSparseImpl(
        IOClass& io,
        IndexType& size, 
        hmemo::HArray<IndexType>& indexes, 
        hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readSparseImpl( size, indexes, static_cast< hmemo::HArray<ValueType>& >( values ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readSparseImpl( io, size, indexes, values );
        }
    }

    static void writeGridImpl( IOClass& io, const hmemo::_HArray& data, const common::Grid& grid )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeGridImpl( static_cast<const hmemo::HArray<ValueType>& >( data ), grid );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeGridImpl( io, data, grid );
        }
    }

    static void readGridImpl( IOClass& io, hmemo::_HArray& data, common::Grid& grid )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readGridImpl( static_cast<hmemo::HArray<ValueType>& >( data ), grid );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readGridImpl( io, data, grid );
        }
    }
};

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
