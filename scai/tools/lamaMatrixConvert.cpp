/**
 * @file tools/lamaMatrixConvert.cpp
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
 * @brief Conversion between matrix file formats, uses FileIO factory
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>

using namespace std;

using namespace scai;
using namespace lama;

static common::ScalarType getType()
{
    common::ScalarType type = common::TypeTraits<double>::stype;

    std::string val;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        scai::common::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );

        if ( env_type == scai::common::ScalarType::UNKNOWN )
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
        cout << "Usage: " << argv[0] << " infile_name outfile_name" << endl;
        cout << "   file format is chosen by suffix, e.g. frm, mtx, txt, psc"  << endl;
        cout << "   --SCAI_TYPE=<data_type> is data type of input file and used for internal representation" << endl;
        cout << "   --SCAI_IO_BINARY=0|1 to force formatted or binary output file" << endl;
        cout << "   --SCAI_IO_TYPE_DATA=<data_type> is data type used for file output" << endl;
        cout << "   " << endl;
        cout << "   Supported types: ";
        vector<common::ScalarType> dataTypes;
        hmemo::_HArray::getCreateValues( dataTypes );

        for ( size_t i = 0; i < dataTypes.size(); ++i )
        {
            cout << dataTypes[i] << " ";
        }

        cout << endl;
        return -1;
    }

    // take double as default

    common::ScalarType type = getType();

    // oops, no factory for storage, only for matrix

    _MatrixPtr matrixPtr( _Matrix::getMatrix( Format::CSR, type ) );
    _Matrix&   matrix = *matrixPtr;

    std::string inFileName = argv[1];

    // use supported file format

    matrix.readFromFile( inFileName );

    cout << "read CSR matrix : " << matrix << endl;

    std::string outFileName = argv[2];

    // use supported file format

    cout << ": write LAMA default format" << endl;

    matrix.writeToFile( outFileName );

    cout << "written matrix to file " << outFileName << endl;
}
