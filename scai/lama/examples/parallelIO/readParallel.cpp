/**
 * @file writeParallel.cpp
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
#include <scai/dmemo.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace std;

using namespace scai;
using namespace lama;
using namespace dmemo;

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

    // oops, no factory for storage, only for matrix

    common::unique_ptr<Matrix> matrixPtr( Matrix::getMatrix( Matrix::CSR, type ) );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    Matrix& matrix = *matrixPtr;

    {
        std::ostringstream inFileName;

        inFileName << comm->getRank() << "." << comm->getSize() << "." << argv[1];

        _MatrixStorage& m = const_cast<_MatrixStorage&>( matrix.getLocalStorage() );

        bool errorFlag = false;

        try
        {
            m.readFromFile( inFileName.str() );
        }
        catch ( common::Exception& e )
        {
            cerr << *comm << ": failed to read " << inFileName.str() << endl;
            errorFlag = true;
        }

        errorFlag = comm->any( errorFlag );

        if ( errorFlag )
        {
            return -1;
        }

        cout << *comm << ": read local part of matrix from file " << inFileName.str() << ": " << m << endl;

        // build the distribution by the sizes

        IndexType globalSize = comm->sum( m.getNumRows() );
        IndexType numColumns = comm->max( m.getNumColumns() );

        // for consistency we have to set the number of columns in each stroage

        dmemo::DistributionPtr rowDist( new dmemo::GenBlockDistribution( globalSize, m.getNumRows(), comm ) );
        dmemo::DistributionPtr colDist( new dmemo::NoDistribution( numColumns ) );

        m.setDimension( m.getNumRows(), numColumns );
        matrix.assign( m, rowDist, colDist );

        cout << *comm << ": distributed matrix = " << matrix << endl;
    }

    dmemo::DistributionPtr dist( new dmemo::CyclicDistribution( matrix.getNumRows(), matrix.getNumRows(), comm ) );

    matrix.redistribute( dist, matrix.getColDistributionPtr() );

    matrix.writeToFile( argv[2] );

    cout << "written CSR matrix : " << matrix << endl;
}
