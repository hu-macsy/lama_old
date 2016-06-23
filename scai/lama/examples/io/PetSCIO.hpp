/**
 * @file PetSCIO.hpp
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
 * @brief Structure that contains IO routines for PETSC
 * @author Thomas Brandes
 * @date 10.06.2016
 */

#pragma once

#include "CRTPFileIO.hpp"

namespace scai
{

namespace lama
{

/** This file format supports the binary format used by PetSC. 
 *
 *   - header information is just at the beginning of the file
 *   - uses CSR format, but the sizes array and not the offsets
 *   - stores data always in BIG endian (x86 has LITTLE endian)
 */

class PetSCIO : 

    public CRTPFileIO<PetSCIO>,         // use type conversions
    public FileIO::Register<PetSCIO>    // register at factory

{

public:

    /** Constructor might reset default values */

    PetSCIO();

    /** Implementation for Printable.:writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    // static method to create an FileIO object for this derived class

    static FileIO* create();

    // registration key for factory

    static std::string createValue();

public:
 
    /** Typed version of writeStorage
     *
     *  This method must be available for implementation of
     *  CRTPFileIO::writeStorage
     */

    template<typename ValueType>
    void writeStorageImpl( const MatrixStorage<ValueType>& storage, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of readStorage */

    template<typename ValueType>
    void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of the writeArray */

    template<typename ValueType>
    void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName )
    __attribute( ( noinline ) );

    /** Typed version of readArray */

    template<typename ValueType>
    void readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName )
    __attribute( ( noinline ) );

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for IO class
};

}

}
