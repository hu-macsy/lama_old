/**
 * @file MetaMatrix.hpp
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
 * SOFTWARE.
 * @endlicense
 *
 * @brief MetaMatrix.hpp
 * @author Kai Buschulte
 * @date 07.05.2012
 * $Id$
 */

#ifndef LAMA_METAMATRIX_HPP_
#define LAMA_METAMATRIX_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/matrix/Matrix.hpp>
#include <lama/distribution/Distribution.hpp>
#include <lama/Context.hpp>
#include <lama/CommunicatorManager.hpp>
#include <lama/CommunicatorFactory.hpp>

// spirit
#include <boost/spirit/include/qi.hpp>

#include <string>
#include <fstream>

namespace lama
{

namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class LAMA_DLL_IMPORTEXPORT MetaMatrix: public Matrix
{
public:
    typedef std::string::const_iterator StringIterator;

    /**
     * @brief Construct a square matrix ... TODO[doxy] Complete Description.
     *
     * @param[in] other    TODO[doxy] Complete Description.
     * @param[in] config   Configuration which will be parsed to generate the defined matrix
     */
    MetaMatrix( Matrix& other, const std::string& config );

    /**
     * @brief virtual destructor of MetaMatrix
     */
    virtual ~MetaMatrix();

    /*
     * @brief If arg is a file it reads the files content. If not it calls the direct interpretation.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     * @param[in] arg     Input string which is a file name or a configuration string
     */
    void interpreteArgument( Matrix& other, const std::string& arg );

    /*
     * @brief Parses the given configuration and constructs the MetaMatrix
     *
     * Parses the given configuration in mConfiguration and constructs the MetaMatrix with the
     * described settings.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     * @param[in] arg     Input string which is a file name or a configuration string
     */
    void parseConfiguration( Matrix& other, const std::string& arg );

    // **** Forward methods to provide matrix functionality ****

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const;

    /* Implementation of pure method of class Matrix. */

    virtual void clear();

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void setIdentity();

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const Matrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& other );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist );

    /* Implementation of pure method of class Matrix. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void getRow( Vector& row, const IndexType globalRowIndex ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void getDiagonal( Vector& diagonal ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Vector& diagonal );

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Scalar scalar );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Vector& scaling );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Scalar scaling );

    /* Implementation of pure method of class Matrix. */

    virtual Scalar getValue( IndexType i, IndexType j ) const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getNumValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void matrixTimesScalar( const Matrix& other, const Scalar alpha );

    /* Implementation of pure method of class Matrix. */

    virtual void matrixTimesMatrix(
        Matrix& result,
        const Scalar alpha,
        const Matrix& B,
        const Scalar beta,
        const Matrix& C ) const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumRows() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumColumns() const;

    /* Implementation of pure method of class Matrix. */

    virtual void setContext( const ContextPtr context );

    /* Implementation of pure method of class Matrix. */

    virtual ContextPtr getContextPtr() const;

    /* Implementation of pure method of class Matrix. */

    virtual MatrixKind getMatrixKind() const;

    /* Implementation of pure method of class Matrix. */

    virtual void prefetch() const;

    /* Implementation of pure method of class Matrix. */

    virtual void wait() const;

    /* Implementation of pure method of class Matrix. */

    virtual void invert( const Matrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual Scalar maxDiffNorm( const Matrix& other ) const;

    /* Implementation of pure method of class Matrix. */

    virtual std::auto_ptr<Matrix> create() const;

    /* Implementation of pure method of class Matrix. */

    virtual std::auto_ptr<Matrix> copy() const;

    /* Implementation of pure method of class Matrix. */

    virtual Scalar::ScalarType getValueType() const;

    /* Implementation of pure method of class Matrix. */

    virtual bool hasDiagonalProperty() const;

    /* Implementation of pure method of class Matrix. */

    virtual void resetDiagonalProperty();

    /* Implementation of pure method of class Matrix. */

    virtual size_t getMemoryUsage() const;

private:
    MatrixPtr mMatrix;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

};

class MatrixConfigGrammar: public qi::grammar<std::string::const_iterator,MatrixPtr( Matrix* ),qi::locals<Matrix*>,
    ascii::space_type>
{
public:
    typedef std::string::const_iterator StringIterator;
    MatrixConfigGrammar();

private:
    qi::rule<StringIterator,MatrixPtr( Matrix* ),qi::locals<Matrix*>,ascii::space_type> mRMatrixConfiguration;

    qi::rule<StringIterator,Matrix*( Matrix* ),ascii::space_type> mRType;

    qi::rule<StringIterator,void( Matrix* ),ascii::space_type> mRContextDef;

    qi::rule<StringIterator,ContextPtr(),ascii::space_type> mRContext;

    qi::rule<StringIterator,ContextManager*(),ascii::space_type> mRContextManager;

    qi::symbols<char,Context::ContextType> mRContextMap;

    qi::rule<StringIterator,void( Matrix* ),ascii::space_type> mRDistributions;

    qi::rule<StringIterator,DistributionPtr( int ),ascii::space_type> mRDistribution;

    qi::rule<StringIterator,CommunicatorPtr(),ascii::space_type> mRComm;

    qi::rule<StringIterator,boost::shared_ptr<CommunicatorManager>(),ascii::space_type> mRCommMan;

    qi::rule<StringIterator,std::string()> mRId;

    MatrixPtr mMatrix;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} //namespace lama

#endif /* LAMA_METAMATRIX_HPP_ */
