/**
 * @file matmul.cpp
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

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Walltime.hpp>

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

    if ( argc != 2 )
    {
        cout << "Usage: " << argv[0] << " infile_name outfile_name" << endl;
        cout << "   file format is chosen by suffix, e.g. frm, mtx, txt, psc"  << endl;
        cout << "   --SCAI_TYPE=<data_type> is data type of input file and used for internal representation" << endl;
        cout << "   --SCAI_IO_BINARY=0|1 to force formatted or binary output file" << endl;
        cout << "   --SCAI_IO_TYPE_DATA=<data_type> is data type used for file output" << endl;
        cout << "   " << endl;
        cout << "   Supported types: ";
        vector<common::scalar::ScalarType> dataTypes;
        hmemo::_HArray::getCreateValues( dataTypes );
        for ( size_t i = 0; i < dataTypes.size(); ++i )
        { 
            cout << dataTypes[i] << " ";
        }
        cout << endl;
        return -1;
    }

    // take double as default 

    common::scalar::ScalarType type = getType();

    cout << "type = " << type << ", ignored, take double" << endl;

    // oops, no factory for storage, only for matrix

    CSRSparseMatrix<double> matrix;

    std::string inFileName = argv[1];

    // use supported file format

    matrix.readFromFile( inFileName );

    cout << "read CSR matrix : " << matrix << endl;

    CSRSparseMatrix<double> cmpMM1( "/home/brandes/MAT/xxt.mat" );

    cout << "cmp mm1 = " << cmpMM1 << endl;

    CSRSparseMatrix<double> cmpMM2( "/home/brandes/MAT/xtx.mat" );

    cout << "cmp mm2 = " << cmpMM2 << endl;

    CSRSparseMatrix<double> matrixT;

    double time = common::Walltime::get();

    matrixT.assignTranspose( matrix );

    time = common::Walltime::get() - time;

    cout << "CSR matrix transposed: " << matrixT << endl;
    cout << "took " << time << " seconds" << endl;

    CSRSparseMatrix<double> mm1;

    time = common::Walltime::get();

    mm1 = matrix * matrixT;

    time = common::Walltime::get() - time;

    cout << "mm1 = " << mm1 << endl;
    cout << "took " << time << " seconds" << endl;

    cmpMM1 -= mm1;

    cout << "Error = " << cmpMM1.maxNorm() << endl;
 
    CSRSparseMatrix<double> mm2;

    time = common::Walltime::get();

    mm2 = matrixT * matrix;

    time = common::Walltime::get() - time;

    cout << "mm2 = " << mm2 << endl;
    cout << "took " << time << " seconds" << endl;

    cmpMM2 -= mm2;

    cout << "Error = " << cmpMM1.maxNorm() << endl;
}
