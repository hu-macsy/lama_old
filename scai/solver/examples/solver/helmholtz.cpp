/**
 * @file solver/examples/solver/helmholtz.cpp
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
 * @brief Example of using double matrices for solvigin complex equation system
 * @author Thomas Brandes
 * @date 30.11.2015
 */

#include "LamaConfig.hpp"
#include "LamaTiming.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Communicator.hpp>
#include <scai/dmemo/distribution/GenBlockDistribution.hpp>

#include <scai/lama/solver/CGNR.hpp>
#include <scai/lama/solver/CG.hpp>
#include <scai/lama/solver/BiCG.hpp>
#include <scai/lama/solver/GMRES.hpp>
#include <scai/lama/solver/TFQMR.hpp>
#include <scai/lama/solver/CGS.hpp>
#include <scai/lama/solver/BiCGstab.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>
#include <scai/lama/solver/criteria/ResidualThreshold.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/lama/norm/L1Norm.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/tracing.hpp>

#include <memory>

/*********************************/

// specify the type of matrix used for solving

typedef ComplexDouble SolverType;

static const double eps = 1e-4;

/*********************************/

using namespace std;

using namespace scai;
using namespace scai::hmemo;
using namespace scai::lama;

/** Count line in a text file */

IndexType countLines( const char* filename )
{
    IndexType lineCount = 0;
    FILE* infile = fopen( filename, "r" );

    if ( infile == NULL )
    {
        COMMON_THROWEXCEPTION( "Could not open file " << filename )
    }

    int ch;

    while ( EOF != ( ch = getc( infile ) ) )
    {
        if ( '\n' == ch )
        {
            ++lineCount;
        }
    }

    cout << "Lines in file " << filename << " : " << lineCount << endl;
    return lineCount;
}

/** Predicate to check whether a type is complex */

template<typename ValueType>
bool isComplexType()
{
    bool isComplex = false;

    switch ( common::TypeTraits<ValueType>::stype )
    {
        case common::ScalarType::DOUBLE_COMPLEX:
        case common::ScalarType::COMPLEX:
        case common::ScalarType::LONG_DOUBLE_COMPLEX :
            isComplex = true;
            break;

        default:
            isComplex = false;
    }

    return isComplex;
}

/** Read in a complex vector of size n.
 *
 *  If ValueType is not complex, 2 * n entries will be read in.
 */

template<typename ValueType>
void readComplexVector( LArray<ValueType>& vec, const IndexType n, const char* fileName )
{
    std::ifstream infile( fileName, ios::in );

    if ( infile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    ContextPtr ctx = Context::getHostPtr();
    bool isComplex = isComplexType<ValueType>();
    // use two times the size if Value type is not complex
    IndexType size = isComplex ? n : 2 * n;
    WriteOnlyAccess<ValueType> wVec( vec, ctx, size );

    for ( IndexType k = 0; k < n; ++k )
    {
        ComplexDouble val;

        if ( infile.fail() )
        {
            COMMON_THROWEXCEPTION( "Could not read values at pos " << k << " in file " << fileName )
        }

        infile >> val;

        if ( isComplex )
        {
            wVec[k] = val;
        }
        else
        {
            wVec[2 * k ]    = val.real();
            wVec[2 * k + 1] = val.imag();
        }

        // cout << "read line " << k << " : " << wVec[2 * k] << ", " << wVec[2 * k + 1 ] << endl;
    }
}

/** Read in a complex matrix of size n x n and nnz non-zero vals in double array of size 2 * n
 *  Each row in input file contains row index, column index, real part, imag part
 */

template<typename ValueType>
void readComplexMatrix( LArray<IndexType>& ia,
                        LArray<IndexType>& ja,
                        LArray<ValueType>& vals,
                        const IndexType n,
                        const IndexType nnz,
                        const char* fileName )
{
    std::ifstream infile( fileName, ios::in );

    if ( infile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    bool isComplex = isComplexType<ValueType>();
    ContextPtr ctx = Context::getHostPtr();
    // use two times the size if Value type is not complex
    IndexType size = isComplex ? nnz : 2 * nnz;
    WriteOnlyAccess<IndexType> wIA( ia, ctx, nnz );
    WriteOnlyAccess<IndexType> wJA( ja, ctx, nnz );
    WriteOnlyAccess<ValueType> wVals( vals, ctx, size );

    for ( IndexType k = 0; k < nnz; ++k )
    {
        if ( infile.fail() )
        {
            COMMON_THROWEXCEPTION( "Could not read values at pos " << k << " in file " << fileName )
        }

        double index;
        ComplexDouble val;
        // be careful: indexes are transfered form 1 : n to 0 : n-1
        infile >> index;
        wIA[k] = static_cast<IndexType>( index - 1 );
        infile >> index;
        wJA[k] = static_cast<IndexType>( index - 1 );
        infile >> val;

        if ( isComplex )
        {
            wVals[ k ] = val;
        }
        else
        {
            wVals[ 2 * k     ] = val.real();
            wVals[ 2 * k + 1 ] = val.imag();
        }

        //      << ", " << wVals[2 * k] << ", " << wVals[2 * k + 1 ] << endl;
        SCAI_ASSERT_LT( wIA[k], n, "illegal row index" )
        SCAI_ASSERT_LT( wJA[k], n, "illegal col index" )
    }
}

/** Conversion of a real matrix to complex matrix
 *  by using the COO format.
 */

template<typename ValueType>
void convertCmplx2Real( LArray<IndexType>& outIA,
                        LArray<IndexType>& outJA,
                        LArray<ValueType>& outValues,
                        const LArray<IndexType>& inIA,
                        const LArray<IndexType>& inJA,
                        const LArray<ValueType>& inValues )
{
    IndexType nnz = inIA.size();
    SCAI_ASSERT_EQUAL( nnz, inJA.size(), "mismatch" )
    SCAI_ASSERT_EQUAL( 2 * nnz, inValues.size(), "mismatch" )
    // size changes from 2 x n to 2 x n, nnz increases by factor 4
    ContextPtr ctx = Context::getHostPtr();
    WriteOnlyAccess<IndexType> wia( outIA, ctx, 4 * nnz );
    WriteOnlyAccess<IndexType> wja( outJA, ctx, 4 * nnz );
    WriteOnlyAccess<ValueType> wvals( outValues, ctx, 4 * nnz );
    ReadAccess<IndexType> ria( inIA, ctx );
    ReadAccess<IndexType> rja( inJA, ctx );
    ReadAccess<ValueType> rvals( inValues, ctx );

    for ( IndexType k = 0; k < nnz; ++k )
    {
        ValueType real = rvals[ 2 * k ];
        ValueType imag = rvals[ 2 * k + 1 ];
        IndexType i = ria[k];
        IndexType j = rja[k];
        wia[4 * k]     = 2 * i;
        wja[4 * k]     = 2 * j;
        wvals[4 * k]   = real;
        wia[4 * k + 1]     = 2 * i;
        wja[4 * k + 1]     = 2 * j + 1;
        wvals[4 * k + 1]   = -imag;
        wia[4 * k + 2]     = 2 * i + 1;
        wja[4 * k + 2]     = 2 * j;
        wvals[4 * k + 2]   = imag;
        wia[4 * k + 3]     = 2 * i + 1;
        wja[4 * k + 3]     = 2 * j + 1;
        wvals[4 * k + 3]   = real;
    }
}

template<typename ValueType>
void testSameVector( const DenseVector<ValueType>& v1, const DenseVector<ValueType>& v2 )
{
    SCAI_ASSERT_EQUAL( v1.size(), v2.size(), "size mismatch" )
    IndexType n = v1.size();
    std::cout << "check " << v1 << " against " << v2 << ", n = " << n << std::endl;
    {
        ContextPtr ctx = Context::getHostPtr();
        ReadAccess<ValueType> r1( v1.getLocalValues(), ctx );
        ReadAccess<ValueType> r2( v2.getLocalValues(), ctx );

        for ( IndexType k = 0; k < n; ++k )
        {
            if ( common::TypeTraits<ValueType>::abs( r1[k] - r2[k] ) > 4.0 *  eps )
            {
                cout << "wrong at pos " << k << ", r1 = " << r1[k] << ", r2 = " << r2[k] << endl;
            }
        }
    }
}

void solveIt( const std::string name, Vector& sol, const _Matrix& mat, const Vector& rhs )
{
    cout << "solveIt, matrix = " << mat << ", rhs = " << rhs << endl;
    LoggerPtr logger( new CommonLogger ( name,
                                         LogLevel::convergenceHistory,
                                         LoggerWriteBehaviour::toConsoleOnly
                                       ) );
    std::unique_ptr<Solver> mySolver( Solver::create( name, name ) );
    mySolver->setLogger( logger );
    IterativeSolver* itSolver = dynamic_cast<IterativeSolver*>( mySolver.get() );
    SCAI_ASSERT( itSolver != NULL, "No iterative solver : " << *mySolver )
    NormPtr norm = NormPtr( new L1Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    CriterionPtr it( new IterationCount( 10000 ) );
    // stop if iteration count reached OR residual threshold is reached
    rt.reset( new Criterion ( it, rt, Criterion::OR ) );
    itSolver->setStoppingCriterion( rt );
    mySolver->initialize( mat );
    double time = common::Walltime::get();
    mySolver->solve( sol, rhs );
    time = common::Walltime::get() - time;
    cout << "Solve took " << time << " seconds." << endl;
}

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        cout << "Call " << argv[0] << " <solver>" << endl;
        return 0;
    }

    // specify the type for which we sole it
    const int n = countLines( "datRHS.txt" );
    const int nnz = countLines( "dataMatrix.txt" );
    cout << "Read in data, n = " << n << ", nnz = " << nnz << endl;
    LArray<SolverType> rhs;
    LArray<SolverType> sol;
    readComplexVector( rhs, n, "datRHS.txt" );
    readComplexVector( sol, n, "datSolution.txt" );
    int size = rhs.size();   // might be 2 * n for double, float
    LArray<IndexType> ia;
    LArray<IndexType> ja;
    LArray<SolverType> mat;
    readComplexMatrix( ia, ja, mat, n, nnz, "dataMatrix.txt" );
    cout << "Read ia : " << ia << endl;
    cout << "Read ja : " << ja << endl;
    cout << "Read vals : " << mat << endl;

    if ( !isComplexType<SolverType>() )
    {
        // convert it to real matrix
        LArray<IndexType> newIA;
        LArray<IndexType> newJA;
        LArray<SolverType> newValues;
        convertCmplx2Real( newIA, newJA, newValues, ia, ja, mat );
        // copy data back
        ia.swap( newIA );
        ja.swap( newJA );
        mat.swap( newValues );
        cout << "new ia : " << ia << endl;
        cout << "new ja : " << ja << endl;
        cout << "new vals : " << mat << endl;
    }

    COOStorage<SolverType> matrix( size, size, ia, ja, mat );
    cout << "COO storage = " << matrix << endl;
    DistributionPtr dist( new NoDistribution( size ) );
    CSRSparseMatrix<SolverType> dMatrix( matrix );  // converts also COO to CSR
    cout << "CSR matrix = " << dMatrix << endl;
    DenseVector<SolverType> rhsVector( rhs, dist );
    DenseVector<SolverType> solVector( sol, dist );
    // just prove that we have got correct input values
    cout << "sol vector = " << solVector << endl;
    DenseVector<SolverType> checkRhsVector = dMatrix * solVector;
    cout << "check rhs = " << checkRhsVector << endl;
    testSameVector( checkRhsVector, rhsVector );
    DenseVector<SolverType> checkSolVector( dist );
    checkSolVector = 1.0;
    ContextPtr ctx = Context::getContextPtr( common::context::Host );
    dMatrix.setContextPtr( ctx );
    checkSolVector.setContextPtr( ctx );
    rhsVector.setContextPtr( ctx );
    dMatrix.setCommunicationKind( _Matrix::SYNCHRONOUS );
    solveIt( argv[1], checkSolVector, dMatrix, rhsVector );
    cout << "Solution available, check for correctness." << endl;
    // testSameVector( checkSolVector, solVector );
}
