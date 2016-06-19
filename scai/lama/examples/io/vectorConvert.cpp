/**
 * @file vectorConvert.cpp
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
 * @brief Conversion between vector file formats
 * @author Thomas Brandes
 * @date 19.06.2016
 */


#include "PetSCIO.hpp"
#include "MatlabIO.hpp"

#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>

using namespace std;

using namespace scai;
using namespace lama;
using namespace hmemo;

static common::scalar::ScalarType getType() 
{
    common::scalar::ScalarType type = common::TypeTraits<double>::stype;
    
    std::string val;
    
    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {   
        scai::common::scalar::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );
        
        if ( env_type == scai::common::scalar::UNKNOWN )
        {   
            std::cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << std::endl;
        }
        
        type = env_type;
    }

    return type;
}

int main( int argc, const char* argv[] )
{
    if ( argc != 3 )
    {
        cout << "Usage: convert infile_name outfile_name" << endl;
        cout << "   file type is chosen by suffix"  << endl;
    }

    // take double as default 

    common::scalar::ScalarType type = getType();

    bool binary = false;   // can be set by environment variable

    common::Settings::getEnvironment( binary, "SCAI_BINARY" );

    PetSCIO petsc_io;

    MatlabIO matlab_io;

    // oops, no factory for storage, only for matrix

    common::unique_ptr<Vector> vectorPtr( Vector::getVector( Vector::DENSE, type ) );

    Vector& vector = *vectorPtr;
    _HArray& array = const_cast<_HArray&>( vector.getLocalValues() );

    std::string inFileName = argv[1];

    if ( _StorageIO::hasSuffix( inFileName, petsc_io.getVectorFileSuffix() ) )
    {
        // use added file format
        petsc_io.readArray( array, inFileName );
    }
    else if ( _StorageIO::hasSuffix( inFileName, matlab_io.getVectorFileSuffix() ) )
    {
        // use added file format
        matlab_io.readArray( array, inFileName );
    }
    else
    {
        // use supported file format
        vector.readFromFile( argv[1] );
    }


    cout << "read array : " << array << endl;

    std::string outFileName = argv[2];

    if ( _StorageIO::hasSuffix( outFileName, petsc_io.getVectorFileSuffix() ) )
    {
        // use added file format
        petsc_io.writeArray( array, outFileName, binary );
    }
    else if ( _StorageIO::hasSuffix( outFileName, matlab_io.getVectorFileSuffix() ) )
    {
        // use added file format
        matlab_io.writeArray( array, outFileName, binary );
    }
    else
    {
        // use supported file format
        vector.writeToFile( outFileName, File::DEFAULT, common::scalar::INTERNAL, binary );
    }
}
