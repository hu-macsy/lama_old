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

#include <scai/lama/io/PartitionIO.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>
#include <scai/dmemo.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

#include "utility.hpp" 

using namespace std;

using namespace scai;
using namespace lama;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        cout << "Usage: " << argv[0] << " infile_name outfile_name [distfile_name]" << endl;
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

    dmemo::DistributionPtr dist;

    if ( argc > 3 )
    {
        dist = PartitionIO::readDistribution( argv[3], comm );
    }

    string inFileName = argv[1];

    bool isPartitioned = false;

    getPartitionFileName( inFileName, isPartitioned, *comm );

    if ( !isPartitioned )
    {
        // read it from a single file

        matrix.readFromFile( inFileName );

        cout << comm << ": read matrix from file " << inFileName << endl;

        if ( dist.get() )
        {
            matrix.redistribute( dist, dist );
        }
    }
    else
    {
        _MatrixStorage& m = const_cast<_MatrixStorage&>( matrix.getLocalStorage() );

        bool errorFlag = false;

        try
        {
            m.readFromFile( inFileName );
        }
        catch ( common::Exception& e )
        {
            cerr << *comm << ": failed to read " << inFileName << endl;
            errorFlag = true;
        }

        errorFlag = comm->any( errorFlag );

        if ( errorFlag )
        {
            return -1;
        }

        cout << *comm << ": read local part of matrix from file " << inFileName << ": " << m << endl;

        if ( dist.get() )
        {
            // we have read a distribution, so it must match

            if ( dist->getLocalSize() != m.getNumRows() )
            {
                errorFlag = true;

                cerr << *comm << ": mismatch distribution and file size" << endl;
            }
        }
        else 
        {
            // we have no distribution so assume a general block distribution

            IndexType globalSize = comm->sum( m.getNumRows() );

            dist.reset( new GenBlockDistribution( globalSize, m.getNumRows(), comm ) );
        }

        // build the distribution by the sizes

        IndexType numColumns = comm->max( m.getNumColumns() );

        // for consistency we have to set the number of columns in each stroage

        m.setDimension( m.getNumRows(), numColumns );

        matrix.assign( m, dist, dist );

        cout << *comm << ": distributed matrix = " << matrix << endl;
    }

    // whatever the distribution may be, we write it in a single file
 
    matrix.writeToFile( argv[2] );

    cout << "written CSR matrix : " << matrix << endl;
}
