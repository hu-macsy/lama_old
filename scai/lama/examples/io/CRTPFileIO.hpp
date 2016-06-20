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

#include <scai/lama/StorageIO.hpp>

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

    void writeStorage( const _MatrixStorage& storage, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage( _MatrixStorage& storage, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::writeArray  */

    virtual void writeArray( const hmemo::_HArray& array, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readArray  */

    virtual void readArray( hmemo::_HArray& array, const std::string& fileName );

    /** Default implementation for removeFile */

    virtual int deleteFile( const std::string& fileName );
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
    static void writeStorageImpl( Derived&, const _MatrixStorage&, const std::string& )
    {
    }

    static void readStorageImpl( Derived&, _MatrixStorage&, const std::string& )
    {
    }

    static void writeArrayImpl( Derived&, const hmemo::_HArray&, const std::string& )
    {
    }

    static void readArrayImpl( Derived&, hmemo::_HArray&, const std::string& )
    {
    }
};

/*
 * Step n
 */
template<class Derived, typename ValueType, typename TailTypes>
struct FileIOWrapper<Derived, common::mepr::TypeList<ValueType, TailTypes> >
{
    static void writeStorageImpl( Derived& io, const _MatrixStorage& storage, const std::string& fileName )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.writeStorageImpl( reinterpret_cast<const MatrixStorage<ValueType>& >( storage ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::writeStorageImpl( io, storage, fileName );
        }
    }

    static void readStorageImpl( Derived& io, _MatrixStorage& storage, const std::string& fileName )
    {
        if ( storage.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readStorageImpl( reinterpret_cast< MatrixStorage<ValueType>& >( storage ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::readStorageImpl( io, storage, fileName );
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
            FileIOWrapper<Derived, TailTypes>::writeArrayImpl( io, array, fileName );
        }
    }

    static void readArrayImpl( Derived& io, hmemo::_HArray& array, const std::string& fileName )
    {
        if ( array.getValueType() == common::getScalarType<ValueType>() )
        {
            io.readArrayImpl( reinterpret_cast< hmemo::HArray<ValueType>& >( array ), fileName );
        }
        else
        {
            FileIOWrapper<Derived, TailTypes>::readArrayImpl( io, array, fileName );
        }
    }
};

/* --------------------------------------------------------------------------------- */
/*    Implementation of methods for CRTPFileIO                                       */
/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeStorage( const _MatrixStorage& storage, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( _StorageIO::hasSuffix( fileName, this->getMatrixFileSuffix() ),
                 fileName << " illegal file name for storage, must have suffix " << getMatrixFileSuffix() )

    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeStorageImpl( ( Derived& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readStorage( _MatrixStorage& storage, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( _StorageIO::hasSuffix( fileName, getMatrixFileSuffix() ),
                 fileName << " illegal, must have suffix " << getMatrixFileSuffix() )

    // just call the corresponding typed routine 

    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::readStorageImpl( ( Derived& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeArray( const hmemo::_HArray& array, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( _StorageIO::hasSuffix( fileName, this->getMatrixFileSuffix() ),
                 fileName << " illegal file name for storage, must have suffix " << getVectorFileSuffix() )

    // now call the corresponding typed routine, use meta-programming to get the correct type
    
    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::writeArrayImpl( ( Derived& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readArray( hmemo::_HArray& array, const std::string& fileName ) 
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( _StorageIO::hasSuffix( fileName, this->getVectorFileSuffix() ),
                 fileName << " illegal file name for array, must have suffix " << this->getVectorFileSuffix() )

    // just call the corresponding typed routine 

    FileIOWrapper<Derived, SCAI_ARITHMETIC_HOST_LIST>::readArrayImpl( ( Derived& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
int CRTPFileIO<Derived>::deleteFile( const std::string& fileName )
{
    int rc = -1;

    if ( _StorageIO::hasSuffix( fileName, this->getMatrixFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else if ( _StorageIO::hasSuffix( fileName, this->getVectorFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else
    {
        SCAI_LOG_WARN( Derived::logger, "do not delete file with unknown suffix" )
    }

    return rc;
}

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
