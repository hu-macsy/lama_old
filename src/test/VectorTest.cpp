/**
 * @file VectorTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class VectorTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 06.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/scoped_array.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/norm/MaxNorm.hpp>
#include <lama/Context.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>

#include <lama/distribution/NoDistribution.hpp>

#include <test/Configuration.hpp>
#include <test/TestSparseMatrices.hpp>
#include <test/EquationHelper.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

#if ARITHMETIC_TYPE_CNT == 1
typedef boost::mpl::list<ARITHMETIC_TYPE0> test_types;
#elif ARITHMETIC_TYPE_CNT == 2
typedef boost::mpl::list<ARITHMETIC_TYPE0, ARITHMETIC_TYPE1> test_types;
#else
typedef boost::mpl::list<ARITHMETIC_TYPE0, ARITHMETIC_TYPE1, ARITHMETIC_TYPE2> test_types;
#endif

/* --------------------------------------------------------------------- */

struct VectorTestConfig
{
    VectorTestConfig()
    {
        m_inputVectorBaseName = Configuration::getInstance().getPath() + "/testVector";
        m_formattedInputVectorBaseName = m_inputVectorBaseName + "Formatted";
        m_xdrDoubleInputVectorBaseName = m_inputVectorBaseName + "XDRDouble";
    }

    ~VectorTestConfig()
    {
    }

    std::string m_inputVectorBaseName;
    std::string m_formattedInputVectorBaseName;
    std::string m_xdrDoubleInputVectorBaseName;
};

BOOST_FIXTURE_TEST_SUITE( VectorTest, VectorTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.VectorTest" )

/* --------------------------------------------------------------------- */

template<typename T>
void verifySameVector( Vector& v1, Vector& v2 )
{
    BOOST_CHECK_EQUAL( v1.getValueType(), v2.getValueType() );

    BOOST_CHECK_EQUAL( v1.size(), v2.size() );

    BOOST_CHECK_EQUAL( v1.getDistribution(), v2.getDistribution() );

    IndexType n = v1.size();

    for ( IndexType i = 0; i < n; ++i )
    {
        // BOOST_CHECK_CLOSE: cannot be used for Complex<T>
        LAMA_CHECK_CLOSE( v1.getValue( i ), v2.getValue( i ), 1 );
    }
}

/* --------------------------------------------------------------------- */

template<typename T>
void verifyVectorWithScalar( Vector& v, Scalar s )
{
    IndexType n = v.size();

    for ( IndexType i = 0; i < n; ++i )
    {
        LAMA_CHECK_CLOSE( v.getValue( i ), s, 1 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, T, test_types )
{
    typedef T ValueType;

    IndexType n = 4;
    DenseVector<ValueType> v( n, 1.0 );
    Scalar scalar = 2.0;

    BOOST_CHECK_EQUAL( n, v.size() );
    for (IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( v.getValue( i ), 1.0 );
    }

    DistributionPtr dist( new NoDistribution( v.size() ) );
    DenseVector<ValueType> v3( v, dist );
    verifySameVector<ValueType>( v, v3 );

    DenseVector<ValueType> v4( v );
    verifySameVector<ValueType>( v4, v );

//Constructor DenseVector( const std::string& filename ) is tested in ReadAndWriteVectorTest
}

/* --------------------------------------------------------------------- */

void cleanupfiles( std::string filename )
{
    std::string prefix = Configuration::getInstance().getPath();
    std::string s = prefix + "/" + filename;
    LAMA_LOG_INFO( logger, "Deleting "<< s << " .vec/.frv ." );
    std::remove( ( s + ".vec" ).c_str() );
    std::remove( ( s + ".frv" ).c_str() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ReadAndWriteVectorTest )
{
    IndexType n = 4;

    DenseVector<double> result( n, 5.0 );

    DenseVector<double> vector( n, 5.0 );
    std::string prefix = Configuration::getInstance().getPath();
    std::string testfilename = "ReadAndWriteVectorTestFile";

    //Write and read FORMATTED
    vector.writeToFile( prefix + "/" + testfilename, File::FORMATTED, File::DOUBLE );

    DenseVector<double> vector2( prefix + "/" + testfilename );
    verifySameVector<double>( vector2, result );
    cleanupfiles( testfilename );

    // write and read BINARY

    std::string fileName = prefix + "/" + testfilename;

    LAMA_LOG_INFO( logger, "write " << vector << " to binary file " << fileName );

    vector.writeToFile( fileName, File::BINARY, File::FLOAT );

    LAMA_LOG_INFO( logger, "Read constructur from binary file " << fileName );

    DenseVector<double> vector3( prefix + "/" + testfilename );
    verifySameVector<double>( vector3, result );
    cleanupfiles( testfilename );

    // write and read XDR
    vector.writeToFile( prefix + "/" + testfilename, File::XDR, File::DOUBLE );

    DenseVector<double> vector5( prefix + "/" + testfilename );
    verifySameVector<double>( vector5, result );
    cleanupfiles( testfilename );

//    TODO: Fix bug in mm_typecode_to_str
//    Write and read MMFile
//    vector.writeToFile( prefix + "/" + testfilename, MATRIX_MARKET, DOUBLE );
//    DenseVector<double> vector6( prefix + "/" + testfilename );
//    verifySameVector<double>( vector6, result );
//    cleanupfiles( testfilename );
}

/* ------------------------------------------------------------------------- */

template<typename mt>
void CtorMatrixExpressionTestmethod()
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    LAMA_LOG_INFO( logger, "CtorMatrixExpressionTestmethod, MatrixType = " << typeid( mt ).name() );

    IndexType n = 4;
    IndexType m = 4;

    DenseVector<ValueType> vectorA( n, 3.0 );
    DenseVector<ValueType> vectorB( n, 2.0 );
    DenseVector<ValueType> vectorC( m, 3.0 );

    Scalar s = 2.0;
    Scalar t = 1.0;

    const CSRSparseMatrix<ValueType> CSRn4m4IdentityMatrix( TestSparseMatrices::n4m4IdentityMatrix<ValueType>() );
    const MatrixType n4m4IdentityMatrix( CSRn4m4IdentityMatrix );

    DenseVector<ValueType> vectorResult( n, 8.0 );
    DenseVector<ValueType> vectorResult2( n, 6.0 );
    DenseVector<ValueType> vectorResult3( n, 3.0 );
    DenseVector<ValueType> vectorResult4( n, 5.0 );
    DenseVector<ValueType> vectorResult5( n, 4.0 );
    DenseVector<ValueType> vectorResult6( n, 1.0 );
    DenseVector<ValueType> vectorResult7( m, 6.0 );
    DenseVector<ValueType> vectorResult8( m, 3.0 );

    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x" );
    DenseVector<ValueType> v1( s * n4m4IdentityMatrix * vectorA );
    verifySameVector<ValueType>( v1, vectorResult2 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x" );
    DenseVector<ValueType> v2( n4m4IdentityMatrix * s * vectorA );
    verifySameVector<ValueType>( v2, vectorResult2 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a" );
    DenseVector<ValueType> v3( n4m4IdentityMatrix * vectorA * s );
    verifySameVector<ValueType>( v3, vectorResult2 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x" );
    DenseVector<ValueType> v4( n4m4IdentityMatrix * vectorA );
    verifySameVector<ValueType>( v4, vectorResult3 );

    /* ************************************************************ */

    LAMA_LOG_INFO( logger, "linear algebra expression: a*x*A" );
    DenseVector<ValueType> v19( s * vectorA * n4m4IdentityMatrix);
    verifySameVector<ValueType>( v19, vectorResult7 );

    LAMA_LOG_INFO( logger, "linear algebra expression: x*A*a" );
    DenseVector<ValueType> v20( vectorA * n4m4IdentityMatrix * s );
    verifySameVector<ValueType>( v20, vectorResult7 );

    LAMA_LOG_INFO( logger, "linear algebra expression: x*a*A" );
    DenseVector<ValueType> v21( vectorA * s * n4m4IdentityMatrix );
    verifySameVector<ValueType>( v21, vectorResult7 );

    LAMA_LOG_INFO( logger, "linear algebra expression: x*A" );
    DenseVector<ValueType> v22( vectorC * n4m4IdentityMatrix );
    verifySameVector<ValueType>( v22, vectorResult8 );

    //AdditionExpressionTest in Constructor

    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x+b*y" );
    DenseVector<ValueType> v5( s * n4m4IdentityMatrix * vectorA + t * vectorB );
    verifySameVector<ValueType>( v5, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x+b*y" );
    DenseVector<ValueType> v6( n4m4IdentityMatrix * s * vectorA + t * vectorB );
    verifySameVector<ValueType>( v6, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a+b*y" );
    DenseVector<ValueType> v7( n4m4IdentityMatrix * vectorA * s + t * vectorB );
    verifySameVector<ValueType>( v7, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x+y" );
    DenseVector<ValueType> v8( s * n4m4IdentityMatrix * vectorA + vectorB );
    verifySameVector<ValueType>( v8, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x+y" );
    DenseVector<ValueType> v9( n4m4IdentityMatrix * s * vectorA + vectorB );
    verifySameVector<ValueType>( v9, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a+y" );
    DenseVector<ValueType> v10( n4m4IdentityMatrix * vectorA * s + vectorB );
    verifySameVector<ValueType>( v10, vectorResult );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x+y" );
    DenseVector<ValueType> v11( n4m4IdentityMatrix * vectorA + vectorB );
    verifySameVector<ValueType>( v11, vectorResult4 );

    //SubstractionExpressionTest in Constructor

    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x-b*y" );
    DenseVector<ValueType> v12( s * n4m4IdentityMatrix * vectorA - t * vectorB );
    verifySameVector<ValueType>( v12, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x-b*y" );
    DenseVector<ValueType> v13( n4m4IdentityMatrix * s * vectorA - t * vectorB );
    verifySameVector<ValueType>( v13, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a-b*y" );
    DenseVector<ValueType> v14( n4m4IdentityMatrix * vectorA * s - t * vectorB );
    verifySameVector<ValueType>( v14, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x-y" );
    DenseVector<ValueType> v15( s * n4m4IdentityMatrix * vectorA - vectorB );
    verifySameVector<ValueType>( v15, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x-y" );
    DenseVector<ValueType> v16( n4m4IdentityMatrix * s * vectorA - vectorB );
    verifySameVector<ValueType>( v16, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a-y" );
    DenseVector<ValueType> v17( n4m4IdentityMatrix * vectorA * s - vectorB );
    verifySameVector<ValueType>( v17, vectorResult5 );

    LAMA_LOG_INFO( logger, "linear algebra expression: A*x-y" );
    DenseVector<ValueType> v18( n4m4IdentityMatrix * vectorA - vectorB );
    verifySameVector<ValueType>( v18, vectorResult6 );

    LAMA_LOG_INFO( logger, "Exceptiontests" );

    DenseVector<ValueType> vec1( 6, 0.0 );

    //Should throw Exception, because of different sizes of matrix and vector
    LAMA_CHECK_THROW( { DenseVector<ValueType> vec( n4m4IdentityMatrix * vec1 ); }, Exception );

}

BOOST_AUTO_TEST_CASE_TEMPLATE( CtorMatrixExpressionTest, T, test_types ) {
    typedef T ValueType;

    CtorMatrixExpressionTestmethod< CSRSparseMatrix<ValueType> >();
    CtorMatrixExpressionTestmethod< ELLSparseMatrix<ValueType> >();
    CtorMatrixExpressionTestmethod< DIASparseMatrix<ValueType> >();
    CtorMatrixExpressionTestmethod< JDSSparseMatrix<ValueType> >();
    CtorMatrixExpressionTestmethod< COOSparseMatrix<ValueType> >();
    CtorMatrixExpressionTestmethod< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CtorVectorExpressionTest, T, test_types )
{
    typedef T ValueType;

    IndexType n = 4;

    DenseVector<ValueType> vectorA( n , 3.0 );
    DenseVector<ValueType> vectorB( n , 5.0 );

    DenseVector<ValueType> vectorResult1( n , 6.0 );
    DenseVector<ValueType> vectorResult2( n , 11.0 );
    DenseVector<ValueType> vectorResult3( n , 26.0 );
    DenseVector<ValueType> vectorResult4( n , 6.0 );
    DenseVector<ValueType> vectorResult5( n , 13.0 );
    DenseVector<ValueType> vectorResult6( n , 8.0 );
    DenseVector<ValueType> vectorResult7( n , -2.0 );
    DenseVector<ValueType> vectorResult8( n , 1.0 );
    DenseVector<ValueType> vectorResult9( n , -7.0 );
    DenseVector<ValueType> vectorResult10( n , -14.0 );

    Scalar s = 2.0;
    Scalar t = 4.0;

    DenseVector<ValueType> v1( s * vectorA );
    verifySameVector<ValueType>( v1, vectorResult1 );

    DenseVector<ValueType> v1a( vectorA * s );
    verifySameVector<ValueType>( v1a, vectorResult1 );

    DenseVector<ValueType> v2 ( s * vectorA + vectorB );
    verifySameVector<ValueType>( v2, vectorResult2 );

    DenseVector<ValueType> v2a ( vectorB + s * vectorA );
    verifySameVector<ValueType>( v2a, vectorResult2 );

    DenseVector<ValueType> v3 ( s * vectorA + t * vectorB );
    verifySameVector<ValueType>( v3, vectorResult3 );

    DenseVector<ValueType> v4( vectorA * s );
    verifySameVector<ValueType>( v4, vectorResult4 );

    DenseVector<ValueType> v5 ( vectorA + s * vectorB );
    verifySameVector<ValueType>( v5, vectorResult5 );

    DenseVector<ValueType> v6 ( vectorA + vectorB );
    verifySameVector<ValueType>( v6, vectorResult6 );

    DenseVector<ValueType> v7 ( vectorA - vectorB );
    verifySameVector<ValueType>( v7, vectorResult7 );

    DenseVector<ValueType> v8 ( s * vectorA - vectorB );
    verifySameVector<ValueType>( v8, vectorResult8 );

    DenseVector<ValueType> v9 ( vectorA - s * vectorB );
    verifySameVector<ValueType>( v9, vectorResult9 );

    DenseVector<ValueType> v10 ( s * vectorA - t * vectorB );
    verifySameVector<ValueType>( v10, vectorResult10 );
}

/* --------------------------------------------------------------------- */

template<typename mt>
void AssignmentOpMatrixExpressionTestmethod( ContextPtr context )
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    IndexType n = 4;
    IndexType m = 4;

    DenseVector<ValueType> vectorA( n, 3.0 );
    DenseVector<ValueType> vectorB( n, 2.0 );
    DenseVector<ValueType> vector( n, 0.0 );
    DenseVector<ValueType> vectorM( m, 0.0 );

    Scalar s = 2.0;
    Scalar t = 1.0;

    const CSRSparseMatrix<ValueType> CSRn4m4IdentityMatrix( TestSparseMatrices::n4m4IdentityMatrix<ValueType>() );
    MatrixType n4m4IdentityMatrix( CSRn4m4IdentityMatrix );

    n4m4IdentityMatrix.setContext( context );

    DenseVector<ValueType> vectorResult( n, 8.0 );
    DenseVector<ValueType> vectorResult2( n, 6.0 );
    DenseVector<ValueType> vectorResult3( n, 3.0 );
    DenseVector<ValueType> vectorResult4( n, 5.0 );
    DenseVector<ValueType> vectorResult5( n, 4.0 );
    DenseVector<ValueType> vectorResult6( n, 1.0 );

    vector = vectorResult;
    verifySameVector<ValueType>( vector, vectorResult );

    Vector& vectorref = vectorResult;

    vector = vectorref;

    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: a*A*x
    LAMA_LOG_INFO( logger, "linear algebra expression: a*A*x" );
    vector = s * n4m4IdentityMatrix * vectorA;
    verifySameVector<ValueType>( vector, vectorResult2 );

    // linear algebra expression: A*a*x
    LAMA_LOG_INFO( logger, "linear algebra expression: A*a*x" );
    vector = n4m4IdentityMatrix * s * vectorA;
    verifySameVector<ValueType>( vector, vectorResult2 );

    // linear algebra expression: A*x*a
    LAMA_LOG_INFO( logger, "linear algebra expression: A*x*a" );
    vector = n4m4IdentityMatrix * vectorA * s;
    verifySameVector<ValueType>( vector, vectorResult2 );

    // linear algebra expression: A*x
    LAMA_LOG_INFO( logger, "linear algebra expression: A*x" );
    vector = n4m4IdentityMatrix * vectorA;
    verifySameVector<ValueType>( vector, vectorResult3 );

    /* ************************************************** */

    // linear algebra expression: a*x*A
    LAMA_LOG_INFO( logger, "linear algebra expression: a*x*A" );
    vectorM = s * vectorA * n4m4IdentityMatrix;
    verifySameVector<ValueType>( vectorM, vectorResult2 );

    // linear algebra expression: x*a*A
    LAMA_LOG_INFO( logger, "linear algebra expression: x*a*A" );
    vectorM = vectorA * s * n4m4IdentityMatrix;
    verifySameVector<ValueType>( vectorM, vectorResult2 );

    // linear algebra expression: x*A*a
    LAMA_LOG_INFO( logger, "linear algebra expression: x*A*a" );
    vectorM = vectorA * n4m4IdentityMatrix * s;
    verifySameVector<ValueType>( vectorM, vectorResult2 );

    // linear algebra expression: x*A
    LAMA_LOG_INFO( logger, "linear algebra expression: x*A" );
    vectorM = vectorA * n4m4IdentityMatrix;
    verifySameVector<ValueType>( vectorM, vectorResult3 );

    //AdditionExpressionTest in Constructor

    // linear algebra expression: a*A*x+b*y
    vector = s * n4m4IdentityMatrix * vectorA + t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: A*a*x+b*y
    vector = n4m4IdentityMatrix * s * vectorA + t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: A*x*a+b*y
    vector = n4m4IdentityMatrix * vectorA * s + t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: a*A*x+y
    vector = s * n4m4IdentityMatrix * vectorA + vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: A*a*x+y
    vector = n4m4IdentityMatrix * s * vectorA + vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: A*x*a+y
    vector = n4m4IdentityMatrix * vectorA * s + vectorB;
    verifySameVector<ValueType>( vector, vectorResult );

    // linear algebra expression: A*x+y
    vector = n4m4IdentityMatrix * vectorA + vectorB;
    verifySameVector<ValueType>( vector, vectorResult4 );

    //SubstractionExpressionTest in Constructor

    // linear algebra expression: a*A*x-b*y
    vector = s * n4m4IdentityMatrix * vectorA - t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: A*a*x-b*y
    vector = n4m4IdentityMatrix * s * vectorA - t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: A*x*a-b*y
    vector = n4m4IdentityMatrix * vectorA * s - t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: a*A*x-y
    vector = s * n4m4IdentityMatrix * vectorA - vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: A*a*x-y
    vector = n4m4IdentityMatrix * s * vectorA - vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: A*x*a-y
    vector = n4m4IdentityMatrix * vectorA * s - vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    // linear algebra expression: A*x-y
    vector = n4m4IdentityMatrix * vectorA - vectorB;
    verifySameVector<ValueType>( vector, vectorResult6 );

    //Exceptiontests

    // Should throw Exception, because of vec1 does not match number of
    // columns in matrix
    // Note: result vector size does not matter as it will be resized

    DenseVector<ValueType> vec1( 6, 1.0 );
    DenseVector<ValueType> vec2( 4, 2.0 );

    LAMA_CHECK_THROW( { vec2 = n4m4IdentityMatrix * vec1; }, Exception );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( AssignmentOpMatrixExpressionTest, T, test_types )
{
    CONTEXTLOOP()
    {
        typedef T ValueType;

        GETCONTEXT( context );

        AssignmentOpMatrixExpressionTestmethod< CSRSparseMatrix<ValueType> >( context );
        AssignmentOpMatrixExpressionTestmethod< ELLSparseMatrix<ValueType> >( context );
        AssignmentOpMatrixExpressionTestmethod< DIASparseMatrix<ValueType> >( context );
        AssignmentOpMatrixExpressionTestmethod< JDSSparseMatrix<ValueType> >( context );
        AssignmentOpMatrixExpressionTestmethod< COOSparseMatrix<ValueType> >( context );
        AssignmentOpMatrixExpressionTestmethod< DenseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

template<typename mt>
void AssignmentVectorExpressionTestmethod( ContextPtr context )
{
    typedef mt ValueType;

    IndexType n = 4;

    DenseVector<ValueType> vectorA( n, 3.0 );
    DenseVector<ValueType> vectorB( n, 5.0 );

    DenseVector<ValueType> vectorResult1( n, 6.0 );
    DenseVector<ValueType> vectorResult2( n, 11.0 );
    DenseVector<ValueType> vectorResult3( n, 26.0 );
    DenseVector<ValueType> vectorResult4( n, 6.0 );
    DenseVector<ValueType> vectorResult5( n, 13.0 );
    DenseVector<ValueType> vectorResult6( n, 8.0 );
    DenseVector<ValueType> vectorResult7( n, -2.0 );
    DenseVector<ValueType> vectorResult8( n, 1.0 );
    DenseVector<ValueType> vectorResult9( n, -7.0 );
    DenseVector<ValueType> vectorResult10( n, -14.0 );

    vectorA.setContext( context );
    vectorB.setContext( context );

    LAMA_LOG_DEBUG( logger, "Using context vecA = " << vectorA.getContext()->getType() );
    LAMA_LOG_DEBUG( logger, "Using context vecB = " << vectorB.getContext()->getType() );

    Scalar s = 2.0;
    Scalar t = 4.0;

    DenseVector<ValueType> vector( n, 0.0 );

    vector = s * vectorA;
    verifySameVector<ValueType>( vector, vectorResult1 );

    vector = vectorA * s;
    verifySameVector<ValueType>( vector, vectorResult1 );

    vector = s * vectorA + vectorB;
    verifySameVector<ValueType>( vector, vectorResult2 );

    vector = vectorB + s * vectorA;
    verifySameVector<ValueType>( vector, vectorResult2 );

    vector = s * vectorA + t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult3 );

    vector = vectorA * s;
    verifySameVector<ValueType>( vector, vectorResult4 );

    vector = vectorA + s * vectorB;
    verifySameVector<ValueType>( vector, vectorResult5 );

    vector = vectorA + vectorB;
    verifySameVector<ValueType>( vector, vectorResult6 );

    vector = vectorA - vectorB;
    verifySameVector<ValueType>( vector, vectorResult7 );

    vector = s * vectorA - vectorB;
    verifySameVector<ValueType>( vector, vectorResult8 );

    vector = vectorA - s * vectorB;
    verifySameVector<ValueType>( vector, vectorResult9 );

    vector = s * vectorA - t * vectorB;
    verifySameVector<ValueType>( vector, vectorResult10 );

    //Exceptiontest

    //Should throw Exception, because of different vector sizes
    DenseVector<ValueType> vec1( 4, 1.0 );
    DenseVector<ValueType> vec2( 6, 2.0 );
    DenseVector<ValueType> vec3( 1, 1.0 );

    LAMA_CHECK_THROW( { vec3 = vec1 + vec2; }, Exception );
    LAMA_CHECK_THROW( { DenseVector<ValueType> vec4( vec1 + 2.0 * vec2 ) ; }, Exception );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( AssignmentVectorExpressionTest, T, test_types )
{
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );

        AssignmentVectorExpressionTestmethod< ValueType >( context );
        AssignmentVectorExpressionTestmethod< ValueType >( context );
        AssignmentVectorExpressionTestmethod< ValueType >( context );
        AssignmentVectorExpressionTestmethod< ValueType >( context );
        AssignmentVectorExpressionTestmethod< ValueType >( context );
        AssignmentVectorExpressionTestmethod< ValueType >( context );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SpecialAssignmentTest, T, test_types )
{
    CONTEXTLOOP()
    {
        typedef T ValueType;

        GETCONTEXT( context );

        IndexType n = 4;

        DenseVector<ValueType> vectorA( n, 3.0 );
        DenseVector<ValueType> vectorB( n, 6.0 );
        DenseVector<ValueType> vectorC( n, 12.0 );

        vectorA.setContext( context );
        vectorB.setContext( context );
        vectorC.setContext( context );

        Scalar s = 2.0;

        LAMA_LOG_INFO( logger, "vectorA *= 2.0" );
        vectorA *= s;
        verifySameVector<ValueType>( vectorA, vectorB );

        vectorA = 6.0;
        LAMA_LOG_INFO( logger, "vectorA += vectorB" );
        vectorA += vectorB;
        verifySameVector<ValueType>( vectorA, vectorC );

        LAMA_LOG_INFO( logger, "vectorA -= vectorB" );
        vectorA -= vectorB;
        verifySameVector<ValueType>( vectorA, vectorB );

        vectorA = 1.0;
        vectorA += 5.0 * vectorA;
        verifySameVector<ValueType>( vectorA, vectorB );

        vectorA = 0.0;
        vectorA += 2.0 * vectorB;
        verifySameVector<ValueType>( vectorA, vectorC );

        DenseVector<ValueType> vectorWrong( n + 1, 6 );

        LAMA_LOG_INFO( logger, "vector(4) += vector(5) should fail" );
        LAMA_CHECK_THROW(
        {   vectorA += vectorWrong;}, Exception );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( operatorDotProductTest, T, test_types )
{
    CONTEXTLOOP()
    {
        typedef T ValueType;

        GETCONTEXT( context );

        IndexType n = 4;

        DenseVector<ValueType> v1( n, 4.0 );
        DenseVector<ValueType> v2( n, 8.0 );

        v1.setContext( context );
        v2.setContext( context );

        Scalar result = v1.dotProduct( v2 );

        BOOST_CHECK_EQUAL( 128.0, result );

        DenseVector<ValueType> v3( n, -2.0 );
        DenseVector<ValueType> v4( n, 14.0 );

        result = v3.dotProduct( v4 );

        BOOST_CHECK_EQUAL( -112.0, result );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MinMaxTest )
{
    DenseVector<double> resVector( 4, 0.0 );

    HostWriteAccess<double> hwares( resVector.getLocalValues() );
    hwares[0] = 9.0;
    hwares[1] = -2.0;
    hwares[2] = 3.0;
    hwares[3] = 6.0;
    hwares.release();

    Scalar s = resVector.min();
    Scalar t = resVector.max();

    BOOST_CHECK_EQUAL( s.getValue<double>(), -2.0 );
    BOOST_CHECK_EQUAL( t.getValue<double>(), 9.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SwapTest )
{
    DenseVector<double> v1( 4, 0.0 );
    DenseVector<double> v2( 4, 1.0 );
    DenseVector<float> v3( 4, 1.0 );

    v1.swap( v2 );

    for ( IndexType i = 0; i < v1.size(); i++ )
    {
        BOOST_CHECK_EQUAL( v1.getValue(i), 1.0 );
        BOOST_CHECK_EQUAL( v2.getValue(i), 0.0 );
    }

    //Should throw exception, because of different vector types
    LAMA_CHECK_THROW( { v1.swap( v3 ); }, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( AssignTest )
{
    DenseVector<float> v1( 4, 0.0 );
    DenseVector<double> v2( 3, 0.0 );
    DenseVector<double> v3( 4, 0.0 );
    DenseVector<float> v4( 4, 1.0 );
    DenseVector<float> v5( 5, 1.0 );

    // Should throw exception, because of different vector sizes

    LAMA_CHECK_THROW( { v1 += v5; }, Exception );

    v3.assign( v4 );
    v1 = v2 = v3;

    BOOST_REQUIRE_EQUAL( v2.size(), v3.size() );

    for ( IndexType i = 0; i < v1.size(); i++ )
    {
        BOOST_CHECK_EQUAL( v1.getValue( i ), 1.0f );
        BOOST_CHECK_EQUAL( v2.getValue( i ), 1.0 );
        BOOST_CHECK_EQUAL( v3.getValue( i ), 1.0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorGetValueTypeTest )
{
    DenseVector<float> v1( 4, 0.0 );
    DenseVector<double> v2( 4, 0.0 );

    BOOST_CHECK_EQUAL( v1.getValueType(), Scalar::FLOAT );
    BOOST_CHECK_EQUAL( v2.getValueType(), Scalar::DOUBLE );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( WriteAtTest )
{
    DenseVector<double> v( 4, 0.0 );
    LAMA_WRITEAT_TEST( v );
}

/* --------------------------------------------------------------------- */

template<typename mt>
void operatorMatrixTimesVectorTestmethod()
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    const IndexType n = 6;
    const IndexType m = 4;

    // Matrix Vector multiplication test 1 ValueType precision
    MatrixType matrixA( TestSparseMatrices::n4m4MatrixA1<ValueType>() );

    DenseVector<ValueType> vectorA( m, 0.0 );
    DenseVector<ValueType> vectorAcalc( matrixA * vectorA );
    DenseVector<ValueType> vectorErg0( m, 0.0 );
    verifySameVector<ValueType>( vectorErg0, vectorAcalc );

    // Matrix Vector multiplication test 2 ValueType precision
    MatrixType matrixA1( TestSparseMatrices::n6m4MatrixD1<ValueType>() );

    DenseVector<ValueType> vectorA1( m, 0.0 );
    DenseVector<ValueType> vectorErg1( n, 0.0 );
    DenseVector<ValueType> vectorA1calc( matrixA1 * vectorA1 );
    verifySameVector<ValueType>( vectorErg1, vectorA1calc );

    // Matrix Vector multiplication test 3 ValueType precision
    MatrixType matrixB1( TestSparseMatrices::n4m6MatrixD2<ValueType>() );

    DenseVector<ValueType> vectorErg2( m, 0.0 );
    DenseVector<ValueType> vectorA2( n, 0.0 );

    DenseVector<ValueType> vectorB1( matrixB1 * vectorA2 );
    verifySameVector<ValueType>( vectorErg2, vectorB1 );

    // Matrix Vector multiplication test 4 ValueType precision
    ValueType vectorBvalues[] =
    { 1.0f, 1.1f, 1.3f, 1.0f };
    DenseVector<ValueType> vectorB( m, vectorBvalues );
    DenseVector<ValueType> vectorBinput( m, 1.0 );

    DenseVector<ValueType> vectorBcalc( matrixA * vectorBinput );
    verifySameVector<ValueType>( vectorB, vectorBcalc );

    // Matrix Vector multiplication test 3 ValueType precision
    DenseVector<ValueType> vectorErg3( m, 0.0 );

    DenseVector<ValueType> vectorA3( m, 0.0 );
    vectorA3 = matrixA * vectorA3;
    verifySameVector<ValueType>( vectorErg3, vectorB1 );

    // Matrix Vector multiplication with non constant input vector
    ValueType vectorAValues[] =
    { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    ValueType vectorErgValues[] =
    { 29.0, 53.0, 36.0, 28.0 };
    DenseVector<ValueType> vectorErg4( m, vectorErgValues );
    DenseVector<ValueType> vectorA4( n, vectorAValues );

    DenseVector<ValueType> vectorB4( matrixB1 * vectorA4 );
    verifySameVector<ValueType>( vectorErg4, vectorB4 );

    MatrixType matrixZ1( TestSparseMatrices::n4m4MatrixA1<ValueType>() );

    //Should throw Exception, because of different sizes of matrix and vector
    // LAMA_CHECK_THROW( { DenseVector<ValueType> d( matrixZ1 * vectorA4 ); }, Exception );

    LAMA_LOG_INFO( logger, "check for exception" );
    DenseVector<ValueType> wrongVectorErg4( m + 1, 1.0 );
    LAMA_CHECK_THROW( { vectorA4 = matrixZ1 * wrongVectorErg4; }, Exception );
    LAMA_LOG_INFO( logger, "check for exception done" );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( operatorMatrixTimeVectorTestold, T, test_types )
{
    typedef T ValueType;

    operatorMatrixTimesVectorTestmethod< CSRSparseMatrix<ValueType> >();
    operatorMatrixTimesVectorTestmethod< ELLSparseMatrix<ValueType> >();
    operatorMatrixTimesVectorTestmethod< DIASparseMatrix<ValueType> >();
    operatorMatrixTimesVectorTestmethod< JDSSparseMatrix<ValueType> >();
    operatorMatrixTimesVectorTestmethod< COOSparseMatrix<ValueType> >();
    operatorMatrixTimesVectorTestmethod< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */

template<typename mt>
void operatorVectorTimesMatrixTestmethod()
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    const IndexType n = 6;
    const IndexType m = 4;

    // Matrix Vector multiplication test 1 ValueType precision
    LAMA_LOG_INFO( logger, "4x4 Matrix" )
    MatrixType matrixA( TestSparseMatrices::n4m4MatrixA1<ValueType>() );

    DenseVector<ValueType> vectorA( m, 0.0 );
    DenseVector<ValueType> vectorAcalc( vectorA * matrixA );
    DenseVector<ValueType> vectorErg0( m, 0.0 );
    verifySameVector<ValueType>( vectorErg0, vectorAcalc );

    // Matrix Vector multiplication test 2 ValueType precision
    LAMA_LOG_INFO( logger, "6x4 Matrix" )
    MatrixType matrixA1( TestSparseMatrices::n6m4MatrixD1<ValueType>() );

    DenseVector<ValueType> vectorA1( n, 0.0 );
    DenseVector<ValueType> vectorErg1( m, 0.0 );
    DenseVector<ValueType> vectorA1calc( vectorA1 * matrixA1 );
    verifySameVector<ValueType>( vectorErg1, vectorA1calc );

    // Matrix Vector multiplication test 3 ValueType precision
    LAMA_LOG_INFO( logger, "4x6 Matrix" )
    MatrixType matrixB1( TestSparseMatrices::n4m6MatrixD2<ValueType>() );

    DenseVector<ValueType> vectorA2( m, 0.0 );
    DenseVector<ValueType> vectorErg2( n, 0.0 );
    DenseVector<ValueType> vectorB1( vectorA2 * matrixB1 );
    verifySameVector<ValueType>( vectorErg2, vectorB1 );

    // Matrix Vector multiplication test 4 ValueType precision
    LAMA_LOG_INFO( logger, "4x4 Matrix" )
    ValueType vectorBvalues[] = { 1.5f, 0.9f, 0.9f, 1.1f };
    DenseVector<ValueType> vectorB( m, vectorBvalues );
    DenseVector<ValueType> vectorBinput( m, 1.0 );

    DenseVector<ValueType> vectorBcalc( vectorBinput * matrixA );
    verifySameVector<ValueType>( vectorB, vectorBcalc );

    // Matrix Vector multiplication test 3 ValueType precision
    LAMA_LOG_INFO( logger, "4x6 Matrix" )
    DenseVector<ValueType> vectorErg3( n, 0.0 );
    DenseVector<ValueType> vectorA3( m, 0.0 );
    vectorA3 = vectorA3 * matrixA;
    verifySameVector<ValueType>( vectorErg3, vectorB1 );

    // Matrix Vector multiplication with non constant input vector
    LAMA_LOG_INFO( logger, "6x4 Matrix" )
    ValueType vectorAValues[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    ValueType vectorErgValues[] = { 38.0, 34.0, 27.0, 45.0 };
    DenseVector<ValueType> vectorErg4( m, vectorErgValues );
    DenseVector<ValueType> vectorA4( n, vectorAValues );

    DenseVector<ValueType> vectorB4( vectorA4 * matrixA1 );
    verifySameVector<ValueType>( vectorErg4, vectorB4 );

    MatrixType matrixZ1( TestSparseMatrices::n4m4MatrixA1<ValueType>() );

    //Should throw Exception, because of different sizes of matrix and vector
    // LAMA_CHECK_THROW( { DenseVector<ValueType> d( matrixZ1 * vectorA4 ); }, Exception );

    LAMA_LOG_INFO( logger, "check for exception" );
    DenseVector<ValueType> wrongVectorErg4( m + 1, 1.0 );
    LAMA_CHECK_THROW( { vectorA4 = matrixZ1 * wrongVectorErg4; }, Exception );
    LAMA_LOG_INFO( logger, "check for exception done" );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( operatorVectorTimesMatrixTestold, T, test_types )
{
    typedef T ValueType;

    operatorVectorTimesMatrixTestmethod< CSRSparseMatrix<ValueType> >();
    operatorVectorTimesMatrixTestmethod< ELLSparseMatrix<ValueType> >();
    operatorVectorTimesMatrixTestmethod< DIASparseMatrix<ValueType> >();
    operatorVectorTimesMatrixTestmethod< JDSSparseMatrix<ValueType> >();
    operatorVectorTimesMatrixTestmethod< COOSparseMatrix<ValueType> >();
    operatorVectorTimesMatrixTestmethod< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( operatorTest, T, test_types )
{
    typedef T ValueType;

    IndexType n = 4;
    DenseVector<ValueType> v( n, 1.0 );

    HostWriteAccess<ValueType> hwa( v.getLocalValues() );

    for (IndexType i = 0; i < n; ++i )
    {
        hwa[i] *= 2.0;
    }

    for (IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( hwa[i], 2.0 );
    }

    hwa.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( xGEMVOperationTest, T, test_types )
{
    typedef T ValueType;
    {
        //y = alpha * A * x + 0 * y

        //nnu = 4 nnc = 3
        const DenseMatrix<ValueType> A( TestSparseMatrices::n4m3MatrixE2<ValueType>() );

        DenseVector<ValueType> x( 3, 1.0 );
        DenseVector<ValueType> resVector( 4, 0.0 );
        DenseVector<ValueType> computeVector( 4, 0.0 );

        HostWriteAccess<ValueType> hwax( x.getLocalValues() );
        hwax[0] = 1.0;
        hwax[1] = 2.0;
        hwax[2] = 3.0;
        hwax.release();

        HostWriteAccess<ValueType> hwares( resVector.getLocalValues() );
        hwares[0] = 9.0;
        hwares[1] = 2.0;
        hwares[2] = 36.0;
        hwares[3] = 6.0;
        hwares.release();

        Scalar alpha = 1.0;

        DenseVector<ValueType> y( 3, 0.0 );

        computeVector = alpha * A * x;

        verifySameVector<ValueType>( resVector, computeVector );

        //constructor test
        DenseVector<ValueType> computeVector2( alpha * A * x );

        verifySameVector<ValueType>( computeVector2, resVector );
    }

    {
        const DenseMatrix<ValueType> A( TestSparseMatrices::n4m4IdentityMatrix<ValueType>() );

        DenseVector<ValueType> x( 4, 1.0 );
        DenseVector<ValueType> y( 4, 1.0 );
        DenseVector<ValueType> vectorRes( 4, 0.0 );

        Scalar scalar = 2.0;
        Scalar zero = 0.0;

        vectorRes = scalar * A * x; // + scalar * y;

        verifyVectorWithScalar<ValueType>( vectorRes, 2.0 );

        //test self-assigment
        x = 2.0;
        x = scalar * A * x;

        verifyVectorWithScalar<ValueType>( x, 4.0 );

        x = 2.0;
        vectorRes = scalar * A * x + scalar * y;

        verifyVectorWithScalar<ValueType>( vectorRes, 6.0 );

        //test constructor
        DenseVector<ValueType> vectorRes2( scalar * A * x + scalar * y );

        verifyVectorWithScalar<ValueType>( vectorRes2, 6.0 );

        //test self-assigment
        x = scalar * A * x + scalar * y;

        verifyVectorWithScalar<ValueType>( x, 6.0 );

        x = 2.0;
        x = scalar * A * x + scalar * x;

        verifyVectorWithScalar<ValueType>( x, 8.0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( xAXPYTest, T, test_types)
{
    typedef T ValueType;

    const DenseMatrix<ValueType> A( TestSparseMatrices::n4m4IdentityMatrix<ValueType>() );

    DenseVector<ValueType> x( 4, 1.0 );
    DenseVector<ValueType> y( 4, 1.0 );
    Scalar scalar = 2.0;

// constructortest
    DenseVector<ValueType> vectorRes( scalar * x + y );
    verifyVectorWithScalar<ValueType>( vectorRes, 3.0 );

    vectorRes = scalar * x + y;
    verifyVectorWithScalar<ValueType>( vectorRes, 3.0 );

// selfassigment test
    x = scalar * x + y;
    verifyVectorWithScalar<ValueType>( x, 3.0 );

    x = 1.0;
    y = scalar * x + y;
    verifyVectorWithScalar<ValueType>( y, 3.0 );

    x = scalar * x + x;
    verifyVectorWithScalar<ValueType>( x, 3.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    DenseVector<double> vector( 4, 2.0 );
    LAMA_WRITEAT_TEST( vector );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
