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
 *      IOWrapper<MatrixMarketIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArray( *this, array );
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
    static void writeStorage( IOClass&, const _MatrixStorage& storage )
    {
        COMMON_THROWEXCEPTION( "writeStorage " << storage << " unsupported, unknown type." )
    }

    static void readStorage( IOClass&, _MatrixStorage& storage )
    {
        COMMON_THROWEXCEPTION( "readStorage " << storage << " unsupported, unknown type." )
    }

    static void writeArray( IOClass&, const hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void writeSparse( IOClass&, const IndexType, const void*, const hmemo::HArray<IndexType>&, const hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void readArray( IOClass&, hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "readArray " << array << " unsupported, unknown type." )
    }

    static void readSparse( IOClass&, IndexType&, void*, hmemo::HArray<IndexType>&, hmemo::_HArray& array )
    {
        COMMON_THROWEXCEPTION( "readSparse unsupported, " << array.getValueType() << " not in SCAI_ARRAY_TYPES"  )
    }

    static void writeGrid( IOClass&, const hmemo::_HArray& data, const common::Grid& )
    {
        COMMON_THROWEXCEPTION( "write " << data << " unsupported, unknown type." )
    }

    static void readGrid( IOClass&, hmemo::_HArray& data, common::Grid& )
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
    static void writeStorage( IOClass& io, const _MatrixStorage& storage )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeStorageImpl( static_cast<const MatrixStorage<ValueType>& >( storage ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeStorage( io, storage );
        }
    }

    static void readStorage( IOClass& io, _MatrixStorage& storage )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readStorageImpl( static_cast< MatrixStorage<ValueType>& >( storage ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readStorage( io, storage );
        }
    }

    static void writeArray( IOClass& io, const hmemo::_HArray& array )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeArrayImpl( static_cast<const hmemo::HArray<ValueType>& >( array ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeArray( io, array );
        }
    }

    static void writeSparse(
        IOClass& io,
        const IndexType size,
        const void* zero,
        const hmemo::HArray<IndexType>& indexes,
        const hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            const ValueType* typedZero = static_cast<const ValueType*>( zero );
            const auto& typedValues = static_cast<const hmemo::HArray<ValueType>& >( values );
            io.writeSparseImpl( size, *typedZero, indexes, typedValues );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeSparse( io, size, zero, indexes, values );
        }
    }

    static void readArray( IOClass& io, hmemo::_HArray& array )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readArrayImpl( static_cast< hmemo::HArray<ValueType>& >( array ) );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readArray( io, array );
        }
    }

    static void readSparse(
        IOClass& io,
        IndexType& size,
        void* zero,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            ValueType* typedZero = static_cast<ValueType*>( zero );
            auto& typedValues = static_cast<hmemo::HArray<ValueType>& >( values );
            io.readSparseImpl( size, *typedZero, indexes, typedValues );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readSparse( io, size, zero, indexes, values );
        }
    }

    static void writeGrid( IOClass& io, const hmemo::_HArray& data, const common::Grid& grid )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeGridImpl( static_cast<const hmemo::HArray<ValueType>& >( data ), grid );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::writeGrid( io, data, grid );
        }
    }

    static void readGrid( IOClass& io, hmemo::_HArray& data, common::Grid& grid )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readGridImpl( static_cast<hmemo::HArray<ValueType>& >( data ), grid );
        }
        else
        {
            IOWrapper<IOClass, TailTypes>::readGrid( io, data, grid );
        }
    }
};

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
