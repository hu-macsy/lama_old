/**
 * @file convert.cpp
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
 * @brief Conversion between matrix file formats, uses FileIO factory
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#include "FileIO.hpp"

#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace std;

using namespace scai;
using namespace lama;

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
    common::Settings::parseArgs( argc, argv );

    if ( argc != 3 )
    {
        cout << "Usage: convert infile_name outfile_name" << endl;
        cout << "   file format is chosen by suffix, e.g. frm, mtx, txt, psc"  << endl;
        cout << "   --SCAI_TYPE=<type> is data type of file input and used for internal data" << endl;
        cout << "   --SCAI_IO_BINARY=0|1 formatted or binary write" << endl;
        cout << "   --SCAI_IO_TYPE_DATA=<type> is data type used for file output" << endl;
        cout << "   " << endl;
        cout << "   Supported data types: float, double, LongDouble, ComplexFloat, ComplexDouble, ComplexLongDouble" << endl;
        exit( -1 );
    }

    // take double as default 

    common::scalar::ScalarType type = getType();

    bool binary = false;   // can be set by environment variable

    common::Settings::getEnvironment( binary, "SCAI_IO_BINARY" );

    // oops, no factory for storage, only for matrix

    common::unique_ptr<Matrix> matrixPtr( Matrix::getMatrix( Matrix::CSR, type ) );

    _MatrixStorage& matrix = const_cast<_MatrixStorage&>( matrixPtr->getLocalStorage() );

    std::string inFileName = argv[1];

    // Note: reading does not care about binary argument, just read as it is

    std::string suffix = FileIO::getSuffix( inFileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        // For reading we expect here that the file data type matches the storage type
        // so SCAI_IO_TYPE_DATA should be ignored for reading 

        fileIO->setDataType( common::scalar::INTERNAL );

        std::cout << "Got from factory: " << *fileIO << std::endl;

        fileIO->readStorage( matrix, inFileName );
    }
    else
    {
        // use supported file format
        matrix.readFromFile( argv[1] );
    }

    cout << "read CSR storage : " << matrix << endl;

    std::string outFileName = argv[2];

    suffix = FileIO::getSuffix( outFileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        std::cout << "Got from factory: " << *fileIO << std::endl;

        // Note. SCAI_IO_TYPE_DATA allows that data is converted

        fileIO->writeStorage( matrix, outFileName );
    }
    else
    {
        // use supported file format
        cout << ": write LAMA default format" << endl;
        matrix.writeToFile( outFileName, File::DEFAULT, common::scalar::INTERNAL, 
                            common::scalar::INDEX_TYPE, common::scalar::INDEX_TYPE, binary );
    }

    cout << "written CSR storage to file " << outFileName << endl;
}