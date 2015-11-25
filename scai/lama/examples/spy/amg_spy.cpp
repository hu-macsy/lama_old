/**
 * @file amg_spy.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Example program to spy sparse structure of a CSR sparse matrices of AMG hierarchy.
 * @author Thomas Brandes
 * @date 17.10.2013
 * @since 1.0.0
 */

#include "Bitmap.hpp"

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/lama/solver/SimpleAMG.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

int main( int argc, char** argv )
{
    CSRSparseMatrix<double> matrix;

    if ( argc < 2 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "spy matrix_filename [ width [ height [ scale ] ] ]" << std::endl;
        exit(1);
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

    int nZoom = 1;

    if ( argc > 4 )
    {
        sscanf( argv[4], "%d",  &nZoom );
    }

    matrix.readFromFile( filename.c_str() );

    // Build matrix hierarchy for Multigrid

    std::string loggerName = "<AMG>";

    LoggerPtr amgLogger( new CommonLogger ( loggerName, scai::lama::LogLevel::completeInformation,
                   LoggerWriteBehaviour::toConsoleOnly,
                   scai::common::shared_ptr<Timer>( new Timer() ) ) );

    scai::common::shared_ptr<SimpleAMG> amgSolver( new SimpleAMG( "SimpleAMG solver", amgLogger ) );

    amgSolver->setHostOnlyLevel( 4 );
    amgSolver->setReplicatedLevel( 5 );
    amgSolver->setMaxLevels( 25 );
    amgSolver->setMinVarsCoarseLevel( 200 );

    amgSolver->initialize( matrix );

    std::cout << "amgSolver has " << amgSolver->getNumLevels() << " levels" << std::endl;

    for ( int level = 0; level < amgSolver->getNumLevels(); ++level )
    {
        const Matrix& mat = amgSolver->getGalerkin( level );
        std::cout << "Galerkin matrix on level " << level << ": " << mat << std::endl;

        HArray<IndexType> ia;
        HArray<IndexType> ja;
        HArray<double> values;

        const _MatrixStorage& local = mat.getLocalStorage();

        local.buildCSRData( ia, ja, values );

        int m = nRows;
        int n = nColumns;
        int s = nZoom;

        if ( m > local.getNumRows() ) 
        {
            m = local.getNumRows();
        }

        if ( n > local.getNumColumns() )
        {
            s *= n / local.getNumColumns();
            n = local.getNumColumns();
        }

        std::cout << "Write png of size " << m << " x " << n << ", zoom = " << s << std::endl;

        Bitmap pic( m, n, s );

        pic.setColor( 240, 120, 0 );  // color for smallest value
        // pic.setColor( 0, 0, 255 );    // color for largest value

        ReadAccess<IndexType> csrIA( ia );
        ReadAccess<IndexType> csrJA( ja );
        ReadAccess<double> csrValues( values );

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

        out_filename << ".galerkin_" << level<< "." << local.getNumRows() << "x" << local.getNumColumns()
                     << "_" << values.size() << ".png";

        pic.write_png_file( out_filename.str().c_str() );

        std::cout << "png files has been written as " << out_filename.str() << std::endl;
    }
 
    for ( int level = 0; level < amgSolver->getNumLevels() - 1; ++level )
    {
        const Matrix& mat = amgSolver->getInterpolation( level );
        std::cout << "Interpolation matrix on level " << level << ": " << mat << std::endl;

        HArray<IndexType> ia;
        HArray<IndexType> ja;
        HArray<double> values;

        const _MatrixStorage& local = mat.getLocalStorage();

        local.buildCSRData( ia, ja, values );

        int m = nRows;
        int n = nRows * local.getNumColumns() / local.getNumRows();
        int s = nZoom;

        if ( m > local.getNumRows() ) 
        {
            m = local.getNumRows();
            n = local.getNumColumns();
            s *= n / local.getNumColumns();
        }

        std::cout << "Write png of size " << m << " x " << n << ", zoom = " << s << std::endl;

        Bitmap pic( m, n, s );

        pic.setColor( 240, 120, 0 );  // color for smallest value
        // pic.setColor( 0, 0, 255 );    // color for largest value

        ReadAccess<IndexType> csrIA( ia );
        ReadAccess<IndexType> csrJA( ja );
        ReadAccess<double> csrValues( values );

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

        out_filename << ".interpolation_" << level<< "." << local.getNumRows() << "x" << local.getNumColumns()
                     << "_" << values.size() << ".png";

        pic.write_png_file( out_filename.str().c_str() );

        std::cout << "png files has been written as " << out_filename.str() << std::endl;
    }
}
