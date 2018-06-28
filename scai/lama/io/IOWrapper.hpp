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

/** Metaprogramming structure to call a routine for each type in a typelist */

template<class Derived, typename TList>
struct IOWrapper;

/*
 * Termination
 */
template<class Derived>
struct IOWrapper<Derived, common::mepr::NullType>
{
    static void writeStorageImpl( Derived&, const _MatrixStorage& storage, const std::string& )
    {
        COMMON_THROWEXCEPTION( "writeStorage " << storage << " unsupported, unknown type." )
    }

    static void readStorageImpl( Derived&, _MatrixStorage& storage, const std::string&, const IndexType, const IndexType )
    {
        COMMON_THROWEXCEPTION( "readStorage " << storage << " unsupported, unknown type." )
    }

    static void writeArrayImpl( Derived&, const hmemo::_HArray& array, const std::string& )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void writeSparseImpl( Derived&, const IndexType, const hmemo::HArray<IndexType>&, const hmemo::_HArray& array, const std::string& )
    {
        COMMON_THROWEXCEPTION( "writeArray " << array << " unsupported, unknown type." )
    }

    static void readArrayImpl( Derived&, hmemo::_HArray& array, const std::string&, const IndexType, const IndexType )
    {
        COMMON_THROWEXCEPTION( "readArray " << array << " unsupported, unknown type." )
    }

    static void readSparseImpl( Derived&, IndexType&, hmemo::HArray<IndexType>&, hmemo::_HArray& array, const std::string& )
    {
        COMMON_THROWEXCEPTION( "readSparse unsupported, " << array.getValueType() << " not in SCAI_ARRAY_TYPES"  )
    }

    static void writeGridImpl( Derived&, const hmemo::_HArray& data, const common::Grid&, const std::string& )
    {
        COMMON_THROWEXCEPTION( "write " << data << " unsupported, unknown type." )
    }

    static void readGridImpl( Derived&, hmemo::_HArray& data, common::Grid&, const std::string& )
    {
        COMMON_THROWEXCEPTION( "read " << data << " unsupported, unknown type." )
    }
};

/*
 * Step n
 */
template<class Derived, typename ValueType, typename TailTypes>
struct IOWrapper<Derived, common::mepr::TypeList<ValueType, TailTypes> >
{
    static void writeStorageImpl( Derived& io, const _MatrixStorage& storage, const std::string& fileName )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeStorageImpl( reinterpret_cast<const MatrixStorage<ValueType>& >( storage ), fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::writeStorageImpl( io, storage, fileName );
        }
    }

    static void readStorageImpl(
        Derived& io,
        _MatrixStorage& storage,
        const std::string& fileName,
        const IndexType offsetRow,
        const IndexType nRows )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readStorageImpl( reinterpret_cast< MatrixStorage<ValueType>& >( storage ), fileName, offsetRow, nRows );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::readStorageImpl( io, storage, fileName, offsetRow, nRows );
        }
    }

    static void writeArrayImpl( Derived& io, const hmemo::_HArray& array, const std::string& fileName )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeArrayImpl( reinterpret_cast<const hmemo::HArray<ValueType>& >( array ), fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::writeArrayImpl( io, array, fileName );
        }
    }

    static void writeSparseImpl( 
        Derived& io, 
        const IndexType size, 
        const hmemo::HArray<IndexType>& indexes, 
        const hmemo::_HArray& values, 
        const std::string& fileName )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeSparseImpl( size, indexes, reinterpret_cast<const hmemo::HArray<ValueType>& >( values ), fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::writeSparseImpl( io, size, indexes, values, fileName );
        }
    }

    static void readArrayImpl(
        Derived& io,
        hmemo::_HArray& array,
        const std::string& fileName,
        const IndexType offset,
        const IndexType n )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readArrayImpl( reinterpret_cast< hmemo::HArray<ValueType>& >( array ), fileName, offset, n );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::readArrayImpl( io, array, fileName, offset, n );
        }
    }

    static void readSparseImpl(
        Derived& io,
        IndexType& size, 
        hmemo::HArray<IndexType>& indexes, 
        hmemo::_HArray& values, 
        const std::string& fileName )
    {
        if ( values.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readSparseImpl( size, indexes, reinterpret_cast< hmemo::HArray<ValueType>& >( values ), fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::readSparseImpl( io, size, indexes, values, fileName );
        }
    }

    static void writeGridImpl( Derived& io, const hmemo::_HArray& data, const common::Grid& grid, const std::string& fileName )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeGridImpl( reinterpret_cast<const hmemo::HArray<ValueType>& >( data ), grid, fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::writeGridImpl( io, data, grid, fileName );
        }
    }

    static void readGridImpl( Derived& io, hmemo::_HArray& data, common::Grid& grid, const std::string& fileName )
    {
        if ( data.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readGridImpl( reinterpret_cast<hmemo::HArray<ValueType>& >( data ), grid, fileName );
        }
        else
        {
            IOWrapper<Derived, TailTypes>::readGridImpl( io, data, grid, fileName );
        }
    }
};

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
