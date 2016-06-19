/**
 * @file CRTPFileIO.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief CRTP wrapper class for all FileIO classes
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#pragma once

#include "FileIO.hpp"

namespace scai
{

namespace lama
{

/* --------------------------------------------------------------------------------- */

/** This help class uses the CRTP pattern to provide for all derived FileIO classes
 *  methods that are implemented in the same way. 
 *
 *  - call typed version of writeStorage and readStorage
 *
 *  @param<tparam> Derived must be a class that is derived from this CRTP class.
 */

template<class Derived>
class CRTPFileIO : public FileIO
{

public:

    /** Implementation of pure virtual method FileIO::writeStorage for all derived classes */

    void writeStorage(
        const _MatrixStorage& storage,
        const std::string& fileName,
        const bool binary = false,
        const common::scalar::ScalarType iaType = common::scalar::INDEX_TYPE,
        const common::scalar::ScalarType jaType = common::scalar::INDEX_TYPE,
        const common::scalar::ScalarType valuesType = common::scalar::INTERNAL ) const;

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage(
        _MatrixStorage& storage,
        const std::string& fileName ) const;

    /** Implementation of pure virtual method FileIO::writeArray  */

    virtual void writeArray(
        const hmemo::_HArray& array,
        const std::string& fileName,
        const bool binary = false,
        const common::scalar::ScalarType valuesType = common::scalar::INTERNAL ) const;

    /** Implementation of pure virtual method FileIO::readArray  */

    virtual void readArray(
        hmemo::_HArray& array,
        const std::string& fileName ) const;
};

/* --------------------------------------------------------------------------------- */
/*  Metaprogramming for different Value Types                                        */
/* --------------------------------------------------------------------------------- */

/** Metaprogramming structure to call a routine for each type in a typelist */

template<class Derived, typename TList>
struct FileIOWrapper;

/*
 * Termination
 */
template<class Derived>
struct FileIOWrapper<Derived, common::mepr::NullType>
{
    static void writeStorageFormatted(
        const _MatrixStorage&,
        const std::string& )
    {
    }

    static void writeStorageBinary(
        const _MatrixStorage&,
        const std::string&,
        const common::scalar::ScalarType,
        const common::scalar::ScalarType,
        const common::scalar::ScalarType )
    {
    }

    static void readStorageTyped(
        _MatrixStorage&,
        const std::string& )
    {
    }

    static void writeArrayFormatted(
        const hmemo::_HArray&,
        const std::string& )
    {
    }

    static void writeArrayBinary(
        const hmemo::_HArray&,
        const std::string&,
        const common::scalar::ScalarType )
    {
    }

    static void readArrayTyped(
        hmemo::_HArray&,
        const std::string& )
    {
    }
};

/*
 * Step n
 */
template<class Derived, typename ValueType, typename TailTypes>
struct FileIOWrapper<Derived, common::mepr::TypeList<ValueType, TailTypes> >
{
    static void writeStorageFormatted(
        const _MatrixStorage& storage,
        const std::string& fileName )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            Derived::writeStorageFormatted( reinterpret_cast<const MatrixStorage<ValueType>& >( storage ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::writeStorageFormatted( storage, fileName );
        }
    }

    static void writeStorageBinary(
        const _MatrixStorage& storage,
        const std::string& fileName,
        const common::scalar::ScalarType iaType,
        const common::scalar::ScalarType jaType,
        const common::scalar::ScalarType valueType )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {   
            Derived::writeStorageBinary( reinterpret_cast<const MatrixStorage<ValueType>& >( storage ), fileName, iaType, jaType, valueType );
        }
        else
        {   
            FileIOWrapper<Derived, TailTypes>::writeStorageBinary( storage, fileName, iaType, jaType, valueType );
        }
    }

    static void readStorageTyped(
        _MatrixStorage& storage,
        const std::string& fileName )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            Derived::readStorageTyped( reinterpret_cast< MatrixStorage<ValueType>& >( storage ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::readStorageTyped( storage, fileName );
        }
    }

    static void writeArrayFormatted(
        const hmemo::_HArray& array,
        const std::string& fileName )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            Derived::writeArrayFormatted( reinterpret_cast<const hmemo::HArray<ValueType>& >( array ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::writeArrayFormatted( array, fileName );
        }
    }

    static void writeArrayBinary(
        const hmemo::_HArray& array,
        const std::string& fileName,
        const common::scalar::ScalarType valueType )
    {   
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {   
            Derived::writeArrayBinary( reinterpret_cast<const hmemo::HArray<ValueType>& >( array ), fileName, valueType );
        }
        else
        {   
            FileIOWrapper<Derived, TailTypes>::writeArrayBinary( array, fileName, valueType );
        }
    }

    static void readArrayTyped(
        hmemo::_HArray& array,
        const std::string& fileName )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            Derived::readArrayTyped( reinterpret_cast< hmemo::HArray<ValueType>& >( array ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::readArrayTyped( array, fileName );
        }
    }
};

/* --------------------------------------------------------------------------------- */
/*    Implementation of methods for CRTPFileIO                                       */
/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeStorage(
    const _MatrixStorage& storage,
    const std::string& fileName,
    const bool binary,
    const common::scalar::ScalarType iaType,
    const common::scalar::ScalarType jaType,
    const common::scalar::ScalarType valueType) const
{
    bool checkedBinary = binary;

    // if binary / formatted is unsupported give a warning and take the other one
    // otherwise the unimplemented versions will throw Exception

    if ( binary )
    {
        if ( ! this->isSupported( binary ) )
        {
            SCAI_LOG_WARN( Derived::logger, "binary unsupported, will write formatted" )
            checkedBinary = false;
        }
    }
    else
    {
        if ( ! this->isSupported( binary ) )
        {
            SCAI_LOG_WARN( Derived::logger, "formatted unsupported, will write binary" )
            checkedBinary = true;
        }
    }

    // now call the corresponding typed routine, use meta-programming to get the correct type

    if ( checkedBinary )
    {
        FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeStorageBinary( storage, fileName, iaType, jaType, valueType );
    }
    else
    {
        FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeStorageFormatted( storage, fileName );
    }
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readStorage(
    _MatrixStorage& storage,
    const std::string& fileName ) const
{
    // just call the corresponding typed routine 

    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::readStorageTyped( storage, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeArray(
    const hmemo::_HArray& array,
    const std::string& fileName,
    const bool binary,
    const common::scalar::ScalarType valueType ) const
{
    bool checkedBinary = binary;
    
    // if binary / formatted is unsupported give a warning and take the other one
    // otherwise the unimplemented versions will throw Exception
    
    if ( binary )
    {   
        if ( ! this->isSupported( binary ) )
        {   
            SCAI_LOG_WARN( Derived::logger, "binary unsupported, will write formatted" )
            checkedBinary = false;
        }
    }
    else
    {   
        if ( ! this->isSupported( binary ) )
        {   
            SCAI_LOG_WARN( Derived::logger, "formatted unsupported, will write binary" )
            checkedBinary = true;
        }
    }
    
    // now call the corresponding typed routine, use meta-programming to get the correct type
    
    if ( checkedBinary )
    {   
        FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeArrayBinary( array, fileName, valueType );
    }
    else
    {   
        FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeArrayFormatted( array, fileName );
    }
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readArray(
    hmemo::_HArray& array,
    const std::string& fileName ) const
{
    // just call the corresponding typed routine 

    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::readArrayTyped( array, fileName );
}

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
