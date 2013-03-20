/**
 * @file MetaMatrix.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * SOFTWARE.using phoenix::val;
 * @endlicense
 *
 * @brief MetaMatrix.cpp
 * @author Kai Buschulte
 * @date 07.05.2012
 * $Id$
 */

// hpp
#include <lama/matrix/MetaMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/SimpleStorageStrategy.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/ContextManager.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/CommunicatorManager.hpp>

// boost
#include <boost/config/warning_disable.hpp>
// spirit
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace lama
{
//        template<typename ValueType>
//        MetaMatrix<ValueType>::logger,
//        "Matrix.MetaMatrix");

LAMA_LOG_DEF_LOGGER( MetaMatrix::logger, "Matrix.MetaMatrix" );

MetaMatrix::MetaMatrix( Matrix& other, const std::string& config )
{
    interpreteArgument( other, config );
}

MetaMatrix::~MetaMatrix()
{
}

void MetaMatrix::interpreteArgument( Matrix& other, const std::string& arg )
{
    std::ifstream configFile;
    configFile.open( arg.c_str() );

    if( configFile )
    {
        LAMA_LOG_DEBUG( logger, "Argument " << arg << " is a file. Reading content now." );
        std::string configuration;

        configFile.unsetf( std::ios::skipws ); // No white space skipping!
        std::copy( std::istream_iterator<char>( configFile ), std::istream_iterator<char>(),
                   std::back_inserter( configuration ) );

        parseConfiguration( other, configuration );
    }
    else
    {
        parseConfiguration( other, arg );
    }
}

void MetaMatrix::parseConfiguration( Matrix& other, const std::string& arg )
{
    LAMA_LOG_INFO( logger, "Parsing configuration " << arg );

    MatrixConfigGrammar configReader; // Our grammar

    using boost::spirit::ascii::space;
    StringIterator first = arg.begin();
    StringIterator last = arg.end();

    bool r = phrase_parse( first, last, configReader( &other ), space, mMatrix );

    if( !r || first != last )
    {
        std::string rest( first, last );
        LAMA_THROWEXCEPTION( "Parsing failure. Stopped at " << rest );
    }

    LAMA_ASSERT( mMatrix, "No root matrix defined in this configuration." );

    LAMA_LOG_DEBUG( logger, "Matrix " << *mMatrix << " is root now." );

}

/***** Grammar Class *****/

LAMA_LOG_DEF_LOGGER( MatrixConfigGrammar::logger, "Matrix.MetaMatrix.Grammar" );

MatrixConfigGrammar::MatrixConfigGrammar()
    : base_type( mRMatrixConfiguration )
{
    using qi::lit;
    using qi::_a;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_r1;
    using qi::_r2;
    using qi::_r3;
    using qi::_val;
    using qi::lexeme;
    using qi::double_;
    using qi::int_;
    using qi::on_error;
    using qi::fail;
    using qi::debug;

    using ascii::char_;

    using phoenix::push_back;
    using phoenix::val;
    using phoenix::construct;
    using phoenix::if_;
    using phoenix::new_;

    mRMatrixConfiguration = mRType( _r1 )[_a = _1] > lit( '{' ) >> -mRContextDef( _a ) >> -mRDistributions( _a )
                            > lit( '}' )[_val = phoenix::construct<MatrixPtr>( _a )];

    mRType = lit( "CSR" )[_val = new_<CSRSparseMatrix<double> >( *_r1 )]
             | lit( "ELL" )[_val = new_<ELLSparseMatrix<double> >( *_r1 )]
             | lit( "JDS" )[_val = new_<JDSSparseMatrix<double> >( *_r1 )]
             | lit( "COO" )[_val = new_<COOSparseMatrix<double> >( *_r1 )]
             | lit( "DIA" )[_val = new_<DIASparseMatrix<double> >( *_r1 )]
             | lit( "Dense" )[_val = new_<DenseMatrix<double> >( *_r1 )] | lit( "SimpleStorageStrategy" )[_val =
                         new_<SimpleStorageStrategy<double> >( *_r1 )];

    /** Context related rules **/

    mRContextDef = ( lit( "context" ) > lit( '=' ) > mRContext > ';' )[phoenix::bind( &Matrix::setContext, _r1, _1 )]
                   | ( lit( "localContext" ) > lit( '=' ) > mRContext > ';' > lit( "haloContext" ) > lit( '=' )
                       > mRContext > ';' )[phoenix::bind( &Matrix::setContext, _r1, _1, _2 )];

    mRContext = mRContextManager[_val = phoenix::bind( &ContextManager::getContext, _1, LAMA_DEFAULT_DEVICE_NUMBER )]
                | ( mRContextManager > '(' > int_ > ')' )[_val = phoenix::bind( &ContextManager::getContext, _1,
                        _2 )];

    mRContextManager = mRContextMap[_val = phoenix::bind( &ContextFactory::getContextManager,
                                           ContextFactory::getFactory(), _1 )];

    mRContextMap.add( "CUDA", Context::CUDA )( "GPU", Context::CUDA )( "Host", Context::Host )( "CPU", Context::Host );

    /** Distribution related rules **/

    mRDistributions =
        ( lit( "rowDistribution=" ) > mRDistribution( phoenix::bind( &Matrix::getNumRows, _r1 ) )
          > lit( ";colDistribution=" )
          > mRDistribution( phoenix::bind( &Matrix::getNumColumns, _r1 ) ) > ';' )[phoenix::bind(
                      &Matrix::redistribute, _r1, _1, _2 )]
        | ( lit( "distribution=" )
            > mRDistribution( phoenix::bind( &Matrix::getNumRows, _r1 ) )[phoenix::bind(
                        &Matrix::redistribute, _r1, _1, _1 )] > lit( ';' ) );

    mRDistribution = ( ( lit( "Block(" ) > mRComm > ")" )[_val = construct<DistributionPtr>(
                           new_<BlockDistribution>( _r1, _1 ) )] | ( lit( "Cyclic(" ) > mRComm > ',' > int_ > ')' )[_val =
                                       construct<DistributionPtr>( new_<CyclicDistribution>( _r1, _2, _1 ) )] );

    int args = 0;
    char** argv = NULL;

    mRComm = qi::eps
             > mRCommMan[_val = phoenix::bind( &CommunicatorManager::getCommunicator,
                                phoenix::bind( &boost::shared_ptr<CommunicatorManager>::get, _1 ),
                                phoenix::ref( args ), phoenix::ref( argv ) )];
    mRCommMan = mRId[_val = phoenix::bind( &CommunicatorFactory::getCommunicatorManager,
                                           CommunicatorFactory::getFactory(), _1 )];

    mRId = char_( "a-zA-Z" )[_val = _1] > *char_( "a-zA-Z" )[_val += _1];

    mRMatrixConfiguration.name( "Configuration" );
    mRComm.name( "Communicator" );
    mRCommMan.name( "CommManager" );
    mRContext.name( "Context" );
    mRContextManager.name( "ContextManager" );
    mRContextDef.name( "Contexts" );
    mRDistribution.name( "Distribution" );
    mRDistributions.name( "Distributions" );
    mRId.name( "ID" );
    mRType.name( "MatrixType" );

    on_error<fail>( mRMatrixConfiguration, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << phoenix::construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << phoenix::construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );
}

/***** Forward Methods MetaMatrix *****/

const char* MetaMatrix::getTypeName() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getTypeName();
}

void MetaMatrix::clear()
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->clear();
}

void MetaMatrix::allocate( const IndexType numRows, const IndexType numColumns )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->allocate( numRows, numColumns );
}

void MetaMatrix::allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->allocate( rowDistribution, colDistribution );
}

void MetaMatrix::setIdentity()
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->setIdentity();
}

void MetaMatrix::assign( const Matrix& other )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->assign( other );
}

void MetaMatrix::assign( const _MatrixStorage& other )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->assign( other );
}

void MetaMatrix::assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->assign( storage, rowDist, colDist );
}

void MetaMatrix::buildLocalStorage( _MatrixStorage& storage ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->buildLocalStorage( storage );
}

void MetaMatrix::redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->redistribute( rowDistribution, colDistribution );
}

void MetaMatrix::getRow( Vector& row, const IndexType globalRowIndex ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->getRow( row, globalRowIndex );
}

void MetaMatrix::getDiagonal( Vector& diagonal ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->getDiagonal( diagonal );
}

void MetaMatrix::setDiagonal( const Vector& diagonal )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->setDiagonal( diagonal );
}

void MetaMatrix::setDiagonal( const Scalar scalar )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->setDiagonal( scalar );
}

void MetaMatrix::scale( const Vector& scaling )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->scale( scaling );
}

void MetaMatrix::scale( const Scalar scaling )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->scale( scaling );
}

Scalar MetaMatrix::getValue( IndexType i, IndexType j ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getValue( i, j );
}

IndexType MetaMatrix::getNumValues() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getNumValues();
}

void MetaMatrix::matrixTimesVector(
    Vector& result,
    const Scalar alpha,
    const Vector& x,
    const Scalar beta,
    const Vector& y ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->matrixTimesVector( result, alpha, x, beta, y );
}

void MetaMatrix::matrixTimesScalar( const Matrix& other, const Scalar alpha )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->matrixTimesScalar( other, alpha );
}

void MetaMatrix::matrixTimesMatrix(
    Matrix& result,
    const Scalar alpha,
    const Matrix& B,
    const Scalar beta,
    const Matrix& C ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->matrixTimesMatrix( result, alpha, B, beta, C );
}

IndexType MetaMatrix::getLocalNumValues() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getLocalNumValues();
}

IndexType MetaMatrix::getLocalNumRows() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getLocalNumRows();
}

IndexType MetaMatrix::getLocalNumColumns() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getLocalNumColumns();
}

void MetaMatrix::setContext( const ContextPtr context )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->setContext( context );
}

ContextPtr MetaMatrix::getContextPtr() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getContextPtr();
}

Matrix::MatrixKind MetaMatrix::getMatrixKind() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getMatrixKind();
}

void MetaMatrix::prefetch() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->prefetch();
}

void MetaMatrix::wait() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->wait();
}

void MetaMatrix::invert( const Matrix& other )
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->invert( other );
}

Scalar MetaMatrix::maxDiffNorm( const Matrix& other ) const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->maxDiffNorm( other );
}

std::auto_ptr<Matrix> MetaMatrix::create() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->create();
}

std::auto_ptr<Matrix> MetaMatrix::copy() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->copy();
}

Scalar::ScalarType MetaMatrix::getValueType() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getValueType();
}

bool MetaMatrix::hasDiagonalProperty() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->hasDiagonalProperty();
}

void MetaMatrix::resetDiagonalProperty()
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    mMatrix->resetDiagonalProperty();
}

size_t MetaMatrix::getMemoryUsage() const
{
    LAMA_ASSERT( mMatrix, "Matrix is not created. Already configured the MetaMatrix?" );
    return mMatrix->getMemoryUsage();
}

} //namespace lama
