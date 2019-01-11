/**
 * @file amg_spy.cpp
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
 * @brief Example program to spy sparse structure of a CSR sparse matrices of AMG hierarchy.
 * @author Thomas Brandes
 * @date 17.10.2013
 */

#include "Bitmap.hpp"

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace scai::solver;

typedef DefaultReal ValueType;

int main( int argc, char** argv )
{
    char image_suffix[] = ".bmp";    // take this as default output format

    // png uses compression but it might happen that it is not supported

    if ( FileIO::canCreate( ".png" ) )
    {
        strcpy( image_suffix, ".png" );
    }

    CSRSparseMatrix<ValueType> matrix;

    if ( argc < 2 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "amg_spy matrix_filename [ width [ height [ scale ] ] ]" << std::endl;
        exit( 1 );
    }

    std::string filename = argv[1];
    int nRows = 800;

    if ( argc > 2 )
    {
        sscanf( argv[2], "%d",  &nRows );
    }

    int nColumns = nRows;

    if ( argc > 3 )
    {
        sscanf( argv[3], "%d",  &nColumns );
    }

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() > 1 )
    {
        std::cout << "This program runs only serial." << std::endl;
        return -1;
    }

    matrix.readFromFile( filename.c_str() );

    // amg solver uses Jacobi solver so we force diagonal property

    matrix.setDiagonalProperty();

    // Build matrix hierarchy for Multigrid
    std::string loggerName = "<AMG>";
    LoggerPtr amgLogger( new CommonLogger ( loggerName, LogLevel::completeInformation,
                                            LoggerWriteBehaviour::toConsoleOnly,
                                            std::shared_ptr<Timer>( new Timer() ) ) );
    std::shared_ptr<SimpleAMG> amgSolver( new SimpleAMG( "SimpleAMG solver", amgLogger ) );
    amgSolver->setHostOnlyLevel( 4 );
    amgSolver->setReplicatedLevel( 5 );
    amgSolver->setMaxLevels( 25 );
    amgSolver->setMinVarsCoarseLevel( 200 );
    amgSolver->initialize( matrix );
    std::cout << "amgSolver has " << amgSolver->getNumLevels() << " levels" << std::endl;

    for ( int level = 0; level < ( int )amgSolver->getNumLevels(); ++level )
    {
        const _Matrix& mat = amgSolver->getGalerkin( level );
        std::cout << "Galerkin matrix on level " << level << ": " << mat << std::endl;
        HArray<IndexType> ia;
        HArray<IndexType> ja;
        HArray<ValueType> values;
        const _MatrixStorage& local = mat.getLocalStorage();
        local.buildCSRData( ia, ja, values );

        IndexType m = nRows;
        IndexType n = nColumns;

        if ( m > local.getNumRows() )
        {
            m = local.getNumRows();
        }

        if ( n > local.getNumColumns() )
        {
            n = local.getNumColumns();
        }

        std::cout << "Write image of size " << m << " x " << n << std::endl;
        Bitmap pic( m, n );
        ReadAccess<IndexType> csrIA( ia );
        ReadAccess<IndexType> csrJA( ja );
        ReadAccess<ValueType> csrValues( values );
        pic.drawCSR( local.getNumRows(), local.getNumColumns(), csrIA.get(), csrJA.get(), csrValues.get() );
        std::ostringstream out_filename;

        if ( filename.find_last_of( "/" ) == std::string::npos )
        {
            out_filename << filename;
        }
        else
        {
            out_filename << filename.substr( filename.find_last_of( "/" ) + 1 );
        }

        out_filename << ".galerkin_" << level << "." << local.getNumRows() << "x" << local.getNumColumns()
                     << "_" << values.size() << image_suffix;
        pic.write( out_filename.str().c_str() );
        std::cout << "image file has been written as " << out_filename.str() << std::endl;
    }

    for ( int level = 0; level < ( int )amgSolver->getNumLevels() - 1; ++level )
    {
        const _Matrix& mat = amgSolver->getInterpolation( level );
        std::cout << "Interpolation matrix on level " << level << ": " << mat << std::endl;
        HArray<IndexType> ia;
        HArray<IndexType> ja;
        HArray<ValueType> values;
        const _MatrixStorage& local = mat.getLocalStorage();
        local.buildCSRData( ia, ja, values );

        IndexType m = nRows;
        IndexType n = nRows * local.getNumColumns() / local.getNumRows();

        if ( m > local.getNumRows() )
        {
            m = local.getNumRows();
            n = local.getNumColumns();
        }

        std::cout << "Write image of size " << m << " x " << n << std::endl;
        Bitmap pic( m, n );
        ReadAccess<IndexType> csrIA( ia );
        ReadAccess<IndexType> csrJA( ja );
        ReadAccess<ValueType> csrValues( values );
        pic.drawCSR( local.getNumRows(), local.getNumColumns(), csrIA.get(), csrJA.get(), csrValues.get() );
        std::ostringstream out_filename;

        if ( filename.find_last_of( "/" ) == std::string::npos )
        {
            out_filename << filename;
        }
        else
        {
            out_filename << filename.substr( filename.find_last_of( "/" ) + 1 );
        }

        out_filename << ".interpolation_" << level << "." << local.getNumRows() << "x" << local.getNumColumns()
                     << "_" << values.size() << image_suffix;
        pic.write( out_filename.str().c_str() );
        std::cout << "image file has been written as " << out_filename.str() << std::endl;
    }
}
