/**
 * @file _Vector.hpp
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
 * @brief Definition of an abstract class for distributed vectors of any kind and any type
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/lama/VectorKind.hpp>

// base classes
#include <scai/dmemo/Distributed.hpp>

// local library
#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/Scalar.hpp>
#include <scai/lama/io/FileIO.hpp>

// others
#include <scai/common/BinaryOp.hpp>
#include <scai/common/CompareOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

#include <scai/common/Factory.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/SCAITypes.hpp>

#include <utility>

namespace scai
{

namespace dmemo
{
class Redistributor;    // forward declaration
}

namespace lama
{

class _Matrix;          // forward declaration

/** Pointer class for a vector, always use of a shared pointer. */

typedef std::shared_ptr<class _Vector> _VectorPtr;

/** Type definition for the key type used for the Vector factory.
 *
 *  The key for vector create is a pair of vector format and the value type.
 */

typedef std::pair<VectorKind, common::ScalarType> VectorCreateKeyType;

/**
 * @brief The class Vector is a abstract type that represents a distributed 1D real or complex vector.
 *
 * Vector is one of the LAMA base types and should be used in all situations where it is not necessary to access a
 * single element or to create a new Vector.
 *
 * As this class is an abstract class, all constructors are protected.
 *
 * The following methods must be implemented by derived classes:
 *
 *  - buildLocalValues to get all local elements on one processor
 *
 * This base class can be used to define dense and sparse vectors of
 * any type.
 */
class COMMON_DLL_IMPORTEXPORT _Vector:

    public common::Factory<VectorCreateKeyType, _Vector*>,
    public dmemo::Distributed

{
public:

    /** @brief More convenient use of the create routine of factory that avoids use of CreateKeyType.
     */
    static _Vector* getVector( const VectorKind format, const common::ScalarType valueType );

    /** @brief More convenient routine to create a dense vector with certain properties.
     *
     *  @param[in] valueType specifies the type of the Vector to be created
     *  @param[in] distribution becomes the distribution of the new vector
     *  @param[in] context optional, becomes the context of the new vector
     */
    static _Vector* getDenseVector(
        const common::ScalarType valueType,
        dmemo::DistributionPtr distribution,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Checks for this vector whether the content of its data is sound.
     *
     * @return false if any of the internal data structures is not okay
     *
     * This method returns the same value on all processors.
     *
     * If any inconsistency has been found an error message should be logged, but it should
     * not throw an exception. This might be done by the caller of this routine to avoid
     * working with inconsistent vectors.
     *
     * \code
     * SCAI_ASSERT_DEBUG( a.isConsistent(), a << ": is invalid matrix after reading" )
     * \endcode
     */
    virtual bool isConsistent() const = 0;

    /**
     * @brief Releases all allocated resources.
     */
    virtual ~_Vector();

    /** Each derived vector must give info about its kind (DENSE or SPARSE). */

    virtual VectorKind getVectorKind() const = 0;

    /**
     * @brief Assigns the values of other to the elements of this.
     *
     * @param[in] other   the vector to get values from.
     * @return            a reference to this.
     */
    _Vector& operator=( const _Vector& other );

    /**
     * @brief Returns the addition of this and other.
     *
     * @param[in] other the vector to do the addition with.
     * @return          a reference to this.
     */
    _Vector& operator+=( const _Vector& other );

    /**
     * @brief Returns the subtraction of this and other.
     *
     * @param[in] other the vector to do the subtraction with.
     * @return          a reference to this.
     */
    _Vector& operator-=( const _Vector& other );

    /**
     * @brief Add a scalar value to all elements of this vector.
     *
     * @param[in] value   the value to add all elements of this with.
     * @return            a reference to this.
     */
    _Vector& operator+=( const Scalar value );

    /**
     * @brief Sub a scalar value to all elements of this vector.
     *
     * @param[in] value   the value to add all elements of this with.
     * @return            a reference to this.
     */
    _Vector& operator-=( const Scalar value );

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value   the value to assign to all elements of this.
     * @return            a reference to this.
     */
    _Vector& operator=( const Scalar value );

    /**
     * @brief Sets the local values of a vector by a dense array.
     *
     * @param[in] values    is the array with all local vector values.
     *
     * The size of the values array must be the same size as the local size of the distribution.
     *
     * Note: Implicit type conversion for the values is supported.
     */

    virtual void setDenseValues( const hmemo::_HArray& values ) = 0;

    /** @brief Allocate and initialize this vector with data from an array
     *
     *  @param[in] values is the array copied to the vector.
     *
     *  Note: the vector is not distributed, i.e. each processor might either set it 
     *        with individual local data or with same data.
     *
     *  \code
     *     HArray<double> arr;
     *     arr.setRandom( 100, 10 ); // fill it with 100 randoms between 0 and 10
     *     Vector& x =
     *     x.setData( arr );         
     *  \endcode
     *
     *  If the vector is redistributed later, it must have been filled with the same values
     *  by each processor.
     */
    void setData( const hmemo::_HArray& values ) 
    {
        allocate( values.size() );
        setDenseValues( values );
    } 

    /** @brief Allocate and initialize this vector with data from an array
     *
     *  @param[in] dist is the distribution of 
     *  @param[in] values become the local values of this vector.
     *
     *  Important: values.size() must be equal to dist->getLocalSize()
     *
     *  \code
     *     IndexType n = 100;
     *     CommunicatorPtr comm = Communicator::getCommunicatorPtr();
     *     DistributionPtr dist( new BlockDistributon( n, comm ) );
     *     HArray<double> arr;     
     *     arr.setRandom( dist->getLocalSize(), 10 ); // every processor fills with random values
     *     DenseVector<double> v;
     *     v.setLocalData( dist, arr );
     *  \endcode
     */
    void setLocalData( dmemo::DistributionPtr dist, const hmemo::_HArray& values ) 
    {
        allocate( dist );
        setDenseValues( values );
    } 

    /** @brief Allocate and initialize this vector with raw values 
     *
     *  @tparam OtherValueType data type of the raw data
     *  @param[in] size becomes the size of the vector and specifies number of entries in values
     *  @param[in] values is pointer to a contiguous array with the raw data
     *
     *  \code
     *    std::vector<float> values;
     *    ....  // build the vector values 
     *    DenseVector<double> v;
     *    v.setRawData( values.size(), &values[0] );
     *  \endcode
     *
     *  Note: the vector is not distributed, i.e. each processor might either set it 
     *        with individual local values or with same values.
     */
    template<typename OtherValueType>
    void setRawData( const IndexType size, const OtherValueType values[] );

    /** Set a replicated vector with sparse vector data
     *
     *  @param[in] n will be the size of the vector
     *  @param[in] nonZeroIndexes positions with non-zero values
     *  @param[in] nonZeroValues values for the non-zero value
     *  @param[in] zeroValue is the 'zero' value, defaults to 0
     *
     *  nonZeroIndexes and nonZeroValues must have the same size. nonZeroIndexes must 
     *  contain valid indexes. They do not have to be sorted.
     */
    void setSparseData( 
        const IndexType n, 
        const hmemo::HArray<IndexType>& nonZeroIndexes, 
        const hmemo::_HArray& nonZeroValues, 
        const Scalar zeroValue = Scalar( 0 ) )
    {
        setSameValue( n, zeroValue );
        fillSparseData( nonZeroIndexes, nonZeroValues, common::BinaryOp::COPY );
    } 

    /** Same as setSparseData but here with raw data for non-zero indexes and values. 
     *
     *  @tparam OtherValueType is the type of the raw data 
     *  @param[in] n will be the size of the vector
     *  @param[in] nnz stands for the number of the non-zero values
     *  @param[in] nonZeroIndexes pointer to array with positions of non-zero values
     *  @param[in] nonZeroValues pointer to array with values
     *  @param[in] zeroValue is the value for all positions that do not appear in nonZeroIndexes
     *
     *  Note: The value type of the raw data might be different to the value type of the vector.
     */
    template<typename OtherValueType>
    void setSparseRawData( 
        const IndexType n, 
        const IndexType nnz,
        const IndexType nonZeroIndexes[],
        const OtherValueType nonZeroValues[],
        const Scalar zeroValue = Scalar( 0 ) );

    /**
     * @brief Sets the local values of a vector by a sparse pattern, i.e. non-zero indexes and values
     *
     * @param[in] nonZeroIndexes   array with all local indexes that have a non-zero entry
     * @param[in] nonZeroValues    array with the values for the nonZeroIndexes
     * @param[in] op               specifies how to deal with available entries, COPY is replace, ADD is sum 
     *
     * Number of non zero indexes and values must be equal, i.e. nonZeroIndexes.size() == nonZeroValues.size()
     *
     * Note: Implicit type conversion for the values is supported. The indexes are local indexes.
     *
     */
    virtual void fillSparseData( 
        const hmemo::HArray<IndexType>& nonZeroIndexes, 
        const hmemo::_HArray& nonZeroValues,
        const common::BinaryOp op ) = 0;

    /**
     * @brief Sets the local data of the vector to zero. 
     *
     * This routine has the same semantic as setValues with an empty array of any type. It is
     * a private routine as it allows a temporary inconsistency between the local part and 
     * the distribution.
     */

    virtual void clearValues() = 0;

    /**
     * This method initilaizes all values of an allocated vector with random numbers.
     *
     * @param[in] bound draw random numbers in the range between 0 and bound (inclusive)
     *
     * For complex vectors a random value is drawn for each real and imaginary part.
     *
     * Keep in mind that bound is an integer value. If you need randonm numbers with other numerical
     * boundaries you should scale them as follows:
     *
     * \code
     *     DistributionPtr dist ( ... );
     *     DenseVector<ValueType> v( dist );
     *     ValueType lb = -1.5, ub = 2.6;
     *     v.fillRandom( 1 );
     *     A = lb + v * ( ub - lb );   // random numbers in the range of lb .. ub
     * \endcode
     */
    virtual void fillRandom( const IndexType bound ) = 0;

    /**
     * This method sets a replicated vector by its size and initializes it with random numbers.
     *
     * In contrary to fillRandom this routine does not require an allocated and maybe uninitialized vector.
     *
     * Be careful: in a parallel environment each processor might initialize the array with different
     * values. By calling Math::srandom( seed ) with the same seed on each processor, it can be forced
     * to have the same values.
     */
    void setRandom( const IndexType n, const IndexType bound )
    {
        allocate ( n );
        fillRandom( bound );
    }

    /**
     * This method sets a distributed vector by its distribution and initializes it with random numbers.
     */
    void setRandom( dmemo::DistributionPtr dist, const IndexType bound )
    {
        allocate ( dist );
        fillRandom( bound );
    }

    /** 
     *  This method gives the vector a size and initializes it with a value.
     *
     *  @param[in] n is the size of the replicated vector
     *  @param[in] value is the value assigned to all elements
     *
     *  \code
     *    DenseVector<double> v1; 
     *    v1.setSameValue( n, value );
     *    DenseVector<double> v2( n );
     *    v2 = value;
     *    DenseVector<double> v3( n, value );
     *  \endcode
     */
    void setSameValue( const IndexType n, const Scalar value )
    {
        allocate( n );
        assign( value );
    }

    /** 
     *  This method gives the vector a distribution and initializes it with a value.
     *
     *  @param[in] dist specifies size of the vector and mapping to the processors
     *  @param[in] value is the value assigned to all elements
     *
     *  \code
     *    DenseVector<double> v1; 
     *    v1.setSameValue( n, value );
     *    DenseVector<double> v2( n );
     *    v2 = value;
     *    DenseVector<double> v3( n, value );
     *  \endcode
     */
    void setSameValue( dmemo::DistributionPtr dist, const Scalar value )
    {
        allocate( dist );
        assign( value );
    }

    /**
     *  Similiar to fillRandom but only replaces the vector elements with a certain probability.
     *
     *  Keep in mind that posititions that are not filled keep their old values. Therefore, in
     *  contrary to fillRandom, the vector must have been initialized before.
     *
     * \code
     *     DistributionPtr dist ( ... );
     *     DenseVector<ValueType> A( dist );
     *     A = 0;
     *     A.fillSparseRandom( 0.5f, 1 );
     * \endcode
     */
    virtual void fillSparseRandom( const float fillRate, const IndexType bound ) = 0;

    /**
     *  Allocate a vector by its size, initialize it with a zero value and fill it sparsely.
     *
     * \code
     *     A.allocate( n );
     *     A = zeroValue;
     *     A.fillSparseRandom( fill, bound );
     * \endcode
     */
    void setSparseRandom( const IndexType n, const Scalar& zeroValue, const float fillRate, const IndexType bound );

    /**
     *  Allocate a vector by its distribution, initialize it with a zero value and fill it sparsely.
     *
     * \code
     *     A.allocate( dist );
     *     A = zeroValue;
     *     A.fillSparseRandom( fill, bound );
     * \endcode
     */
    void setSparseRandom( dmemo::DistributionPtr dist, const Scalar& zeroValue, const float fillRate, const IndexType bound );

    /** @brief This method scales all vector values with a scalar.
     *
     * @param[in] scaling   is the source value.
     *
     * \deprecated{Please use operator*= or operator/= of Vector<ValueType>}
     */
    virtual void _scale( const Scalar scaling );

    /**
     * @brief Elementwise multiplication with another vector (same size), i.e. this[i] = this[i] * other[i]
     *
     * @param[in] other   the vector to multiply to do the multiplication per element
     *
     * \deprecated{Please use Vector<ValueType>::cwiseProdouct}
     */
    virtual void _cwiseProduct( const _Vector& other );

    /**
     * @brief Elementwise disision with another vector (same size), i.e. this[i] = this[i] / other[i]
     *
     * @param[in] other   the vector used for elementwise division
     *
     * \deprecated{Please use Vector<ValueType>::cwiseDivision}
     */
    virtual void _cwiseDivision( const _Vector& other );

    /**
     * This method sets a vector by reading its values from one or multiple files.
     *
     * @param[in] fileName      the filename to read from
     * @param[in] distribution  optional, if set it is the distribution of the vector
     *
     *   \code
     *      DenseVector<double> vector;
     *      vector.readFromFile( "vector.mtx" )                    ! vector only on processor 0
     *      vector.readFromFile( "vector_%r.mtx" )                 ! general block distributed vector, each processor reads it own file
     *      vector.readFromFile( "vector.mtx", rowDist )           ! each processor gets its local part of the vector in one file
     *      vector.readFromFile( "vector_%r.mtx", rowDist )        ! read a partitioned vector with the given distribution
     *   \endcode
     */
    void readFromFile( const std::string& fileName, dmemo::DistributionPtr distribution = dmemo::DistributionPtr() );

    /**
     *  This method sets a vector a reading its values from one or multiple files and also the distribution from a file
     *
     * @param[in] vectorFileName the single or partitioned filename to read from
     * @param[in] distributionFileName the single or partitioned filename with the row distribution of the vector
     *
     *   \code
     *      CSRSparseMatrix<double> vector;
     *      vector.readFromFile( "vector.mtx", "owners.mtx" )
     *      vector.readFromFile( "vector_%r.mtx", "owners.mtx" )
     *      vector.readFromFile( "vector.mtx", "rows%r.mtx" )
     *      vector.readFromFile( "vector_%r.mtx", "rows%r.mtx" )
     *   \endcode
     */
    void readFromFile( const std::string& vectorFileName, const std::string& distributionFileName );

    /**
     * @brief write the vector to an output file
     *
     * @param[in] fileName is the name of the output file (suffix must be added according to the file type)
     * @param[in] fileType format of the output file ("frv" for SAMG, "mtx" for _MatrixMarket), default is to decide by suffix
     * @param[in] dataType representation type for output values, if set it overrides IO settings
     * @param[in] fileMode can be BINARY or FORMATTED, DEFAULT_MODE keeps default/environment settings
     */
    void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::ScalarType dataType = common::ScalarType::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    /**
     * @brief Queries the value type of the vector elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::ScalarType getValueType() const = 0;

    /**
     * @brief Concatenate multiple vectors to a new vector.
     *
     * @param[in] dist specifies the distribution of the concatenated vector.
     * @param[in] vectors is a vector with const pointers/references to the concatenated vectors
     *
     * Note: dist.getGlobalSize() == v[0]->size() + ... v[n-1]->size() 
     *
     * This routine should also be able to deal with aliases, i.e. one ore more of the pointers might be
     * this vector itself.
     */
    virtual void concatenate( dmemo::DistributionPtr dist, const std::vector<const _Vector*>& vectors ) = 0;

    /**
     * @brief Concatenate two vectors to a new vector.
     *
     * @param[in] v1 first part of the new vector
     * @param[in] v2 second part of the new vector
     */
    virtual void cat( const _Vector& v1, const _Vector& v2 );

    /**
     *  @brief Creates a new Vector of the same kind and value type, and same context
     *
     *  /code
     *    const Vector& old = ...
     *    ....
     *
     *  The new vector is a zero vector, neither allocated, nor initialized.
     */
    virtual _Vector* newVector() const = 0;

    /**
     *  @brief copy is a virtual call of the copy constructor of the derived classes
     *
     *  /code
     *    const Vector& old = ...
     *    Vector* new = old.cooy()
     *
     *    // is same as
     *
     *    Vector* new = Vector::create( old.getCreateValue() );
     *    *new = old;
     *  /endcode
     */
    virtual _Vector* copy() const = 0;

    /**
     * @brief Returns the size of the vector.
     *
     * @return  the size of this vector.
     */
    inline IndexType size() const;

    /**
     * @brief Swaps the content of this vector with another vector.
     *
     * @param[in,out] other   the Vector to swap the contents with.
     *
     * Swap is only possible if both vectors are of the same format (DENSE) and
     * have the same value type.
     */
    virtual void swap( _Vector& other ) = 0;

    /** Override default implementation of Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /**
     *  @brief Assigns an arbitrary vector to this vector.
     *
     *  @param[in] other is any vector
     *
     *  This vector inherits the size / distribution from the other vector, but not the context.
     *  This operation will support all kind of conversions, i.e. type (float, double, .. )
     *  and kind ( SPARSE, DENSE ) conversions.
     *
     *  Each vector class has to implement this pure method.
     */
    virtual void assign( const _Vector& other ) = 0;

    /**
     *  Assignment to vector by local values and distribution.
     */
    void assign( const hmemo::_HArray& localValues, dmemo::DistributionPtr distribution );

    /**
     *  Define a non-distributed vector by an array with all its values.
     *
     *  Note: for a correct replication all processors must set the same values.
     */
    void assign( const hmemo::_HArray& globalValues );

    /**
     *  Build an array with all local values of a distributed vector.
     *
     *  @param[in,out] localValues   will be an array that contains local values of the vector
     *  @param[in]     op            specifies how to combine with existing values in localValues
     *  @param[in]     prefLoc       is the location where the values are needed
     *
     *  For different value types, implicit format conversion will be done.
     *  A sparse vector might generate an array with all local values.
     *
     *  If op is not COPY, the binary operation op is applied to existing values in localValues. In this
     *  case, the array localValues must have the local size of the distribution.
     */
    virtual void buildLocalValues( 
        hmemo::_HArray& localValues, 
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const = 0;

    /**
     *  Gather certain local values of a distributed vector.
     *
     *  @param[in,out] localValues   array for the gather values
     *  @param[in]     localIndexes  are the indexes to be gathered
     *  @param[in]     op            specifies how to combine with existing values in localValues
     *  @param[in]     prefLoc       is the location where the values are needed
     *
     *  For different value types, implicit format conversion will be done.
     *
     *  If op is not COPY, the binary operation op is applied to existing values in localValues. In this
     *  case, the array localValues must have the same size as localIndexes.
     */
    virtual void gatherLocalValues( 
        hmemo::_HArray& localValues, 
        const hmemo::HArray<IndexType>& localIndexes,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const = 0;

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value   the value to assign to all elements of this.
     */
    virtual void assign( const Scalar value ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX + scalarBeta * vectorY 
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorPlusVector( const Scalar& alphaS, const _Vector& x, const Scalar& betaS, const _Vector& y ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX * vectorY
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorTimesVector( const Scalar& alphaS, const _Vector& x, const _Vector& y ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX * scalarBeta
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorPlusScalar( const Scalar& alphaS, const _Vector& x, const Scalar& betaS ) = 0;

    /**
     *  @brief Update this vector with another vector elementwise
     * 
     *  @param[in] other is the input vector for setting, must have same distribution
     *  @param[in] op specifies the binary operation for the update
     *  @param[in] swapArgs if true the arguments of the binary operator are swapped
     * 
     *  The call v1.setVector( v2, op ) is equivalent to the following code:
     *
     *  \code
     *      SCAI_ASSERT_EQ_ERROR(( v1.getDistribution(), v2.getDistribuiton(), "mismatch" )
     *      for ( IndexType i = 0; i < v1.size(); ++i )
     *      {
     *          v1[i] = v1[i] op v2[i];    // swapArgs = false
     *          v1[i] = v2[i] op v1[i];    // swapArgs = true
     *      }
     *  \endcode
     *
     *  In contrary to the loop, it can be assumed that the vector operation is full parallel.
     */
    virtual void setVector( const _Vector& other, common::BinaryOp op, const bool swapArgs = false ) = 0;

    /**
     *  @brief Update this vector with a scalar value elementswise
     * 
     *  @param[in] value is th scalar element used for the operation
     *  @param[in] op specifies the binary operation for the update
     *  @param[in] swapScalar if true the operands are swapped
     * 
     *  The call v.setScalar( s, op ) is equivalent to the following code:
     *
     *  \code
     *      for ( IndexType i = 0; i < v1.size(); ++i )
     *      {
     *          v[i] = v[i] op s;    // swapScalar = false
     *          v[i] = s op v[i];    // swapScalar = true
     *      }
     *  \endcode
     *
     *  Here are some examples how this method is used:
     *  \code
     *      v.invert()      v.setScalar( Scalar( 1 ), BinaryOp::DIVIDE, true );
     *      v += s;         v.setScalar( s, BinaryOp::ADD, false );
     *      v *= s;         v.setScalar( s, BinaryOp::MULT, false );
     *  \endcode
     */
    virtual void setScalar( const Scalar value, common::BinaryOp op, const bool swapScalar = false ) = 0;

    /**
     *  @brief Apply a UnaryOp operation for each element of the vector.
     */
    virtual void applyUnary( common::UnaryOp op ) = 0;

    /**
     *  @brief Boolean reduction returns true if all elements fullfill the compare operation with a scalar.
     */
    virtual bool all( common::CompareOp op, const Scalar value ) const = 0;

    /**
     *  @brief Boolean reduction returns true if elementwise comparison with other vector is true for all elements
     */
    virtual bool all( common::CompareOp op, const _Vector& other ) const = 0;

    /**
     * @brief Starts a prefetch to make this valid at the passed context.
     *
     * @param[in] context specifies the location to make this vector valid at
     */
    virtual void prefetch( const hmemo::ContextPtr context ) const = 0;

    /**
     * @brief Starts a prefetch to make data valid at the context of the vector.
     *
     */

    void prefetch() const;

    /**
     * @brief Waits for a possibly running prefetch.
     */
    virtual void wait() const = 0;

    /**
     * @brief Sets the 'preferred' context where data resides and computations are done.
     */
    void setContextPtr( hmemo::ContextPtr context );

    /**
     * @brief Getter function for the context (pointer) of a vector.
     */
    inline hmemo::ContextPtr getContextPtr() const;

    /**
     * @brief Getter function for the context (const reference) of a vector.
     */
    inline const hmemo::Context& getContext() const;

    /**
     * @brief Returns the global memory that is allocated to hold this vector.
     * For a distributed vector all partitions are summed together.
     *
     * @return the memory consumption of this vector.
     */
    virtual size_t getMemoryUsage() const = 0;

    /**
     *  @brief Allocates or reallocates this vector for a given distribution.
     *
     *  All elements of the vector are undefined after this operation.
     *  This operation will allocate memory at the context of this vector.
     */
    virtual void allocate( dmemo::DistributionPtr distributionPtr ) = 0;

    /**
     *  @brief Set distribution and context in one method for convenience
     *
     *  Note: vector itself remains uninitialized.
     */
    void setSpace( dmemo::DistributionPtr distributionPtr, hmemo::ContextPtr ctx )
    {
        // set the context before the allocate so it will be allocated at the right place
        setContextPtr( ctx );
        allocate( distributionPtr );
    }

    /**
     *  @brief Allocates or reallocates this vector as replicted with the given size.
     *
     *  All elements of the vector are undefined after this operation.
     */
    virtual void allocate( const IndexType n ) = 0;

    /**
     * @brief Redistributes this vector to the new passed distribution.
     *
     * @param[in] distribution   the new distribution for this vector.
     *
     * The global vector itself remains unchanged; only local parts
     * can be different now.
     */
    virtual void redistribute( dmemo::DistributionPtr distribution ) = 0;

    /**
     *  @brief Redistribute this vector with a redistributor 
     *
     *  Note: redistributor.getSourceDistribution() == this->getDistribution(),
     *        must hold before the call, this->getDistribution() == redistributor.getTargetDistribution()
     *        is valid after the call.
     */
    virtual void redistribute( const dmemo::Redistributor& redistributor ) = 0;

    /** 
     * @brief Replicate this vector, i.e. redistribute with NoDistribution( size() )
     */
    void replicate();

    /**
     * @brief This method inverts all elements of the vector and is completely local.
     */
    void invert();

    /**
     *  Build the conjugate vector in place.
     */
    void conj();

    /**
     *  Build the absolute in place.
     */
    void abs();

    /**
     *  Calculates the exponentional function of the vector elements in place.
     */
    void exp();

    /**
     *  Calculates the logarithm of the vector elements in place.
     */
    void log();

    /**
     *  Calculates the floor function of the vector elements in place.
     */
    void floor();

    /**
     *  Calculates the ceil function of the vector elements in place.
     */
    void ceil();

    /**
     *  Calculates the square root of the vector elements.
     */
    void sqrt();

    /**
     *  Calculates the sinus of the vector elements.
     */
    void sin();

    /**
     *  Calculates the cosinus of the vector elements.
     */
    void cos();

    /**
     *  Calculates the tangens of the vector elements.
     */
    void tan();

    /**
     *  Calculates the arcus tangens of the vector elements.
     */
    void atan();

    /**
     *  Calculates the pow function for the vector elements with the elements of another vector.
     */
    void powBase( const _Vector& other );

    /**
     *  Calculates the pow function for the vector elements with the elements of another vector.
     */
    void powExp( const _Vector& other );

    /**
     *  Calculates the pow function for a base the vector elements as exponents.
     */
    void powBase( const Scalar base );

    /**
     *  Calculates the pow function for the vector elements as base and an exponent.
     */
    void powExp( const Scalar exp );

protected:

    /**
     *  Constructor of _Vector for derived classes by size and/or context
     */
    explicit _Vector( const IndexType size = 0, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context       is optional, will be Host context.
     */
    explicit _Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Creates a copy of the passed _Vector.
     *
     * @param[in] other   the Vector to take a copy from.
     *
     * Inherits size/distribution, context and content of the passed vector
     */
    _Vector( const _Vector& other );

    /**
     *  @brief Swaps member variables of Vector class.
     */
    void swapVector( _Vector& other );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    hmemo::ContextPtr mContext; //!< decides about location of vector operations

    /** write only the local data to a file, no communication here */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileIO::FileMode fileMode ) const = 0;

    /** write the whole vector into a single file, can imply redistribution */

    void writeToSingleFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** same as writeLocalToFile but also communication for error handling */

    void writeToPartitionedFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    void readFromSingleFile( const std::string& fileName );

    /** Read only the local part from a file, no communication here.
     *
     *  This routine is implemented individually by sparse and dense vectors.
     *
     *  @param[in] fileName is the name of the input file containing the local vector data
     *  @param[in] first is index of first element to read
     *  @param[in] size number of elements to read, if nIndex read up to maximal size
     *  @return    the size of the local vector read in
     *
     *  This routine is private as it allows a temporary inconsistency between the size of 
     *  the local vector data and the distribution.
     */
    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = nIndex ) = 0;

    /** In this version each processor reads from input file its local part. */

    void readFromSingleFile( const std::string& fileName, dmemo::DistributionPtr dist );

    void readFromPartitionedFile( const std::string& myPartitionFileName, dmemo::DistributionPtr dist );
};

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

IndexType _Vector::size() const
{
    return getDistribution().getGlobalSize();
}

hmemo::ContextPtr _Vector::getContextPtr() const
{
    return mContext;
}

const hmemo::Context& _Vector::getContext() const
{
    // Note: mContext is NEVER zero pointer
    return *mContext;
}

template<typename OtherValueType>
void _Vector::setRawData( const IndexType size, const OtherValueType values[] )
{
    allocate( size );

    // use heterogeneous array reference to avoid copy of the raw data

    hmemo::HArrayRef<OtherValueType> valuesArrayRef( size, values );
    setDenseValues( valuesArrayRef );
}

template<typename OtherValueType>
void _Vector::setSparseRawData( 
    const IndexType n, 
    const IndexType nnz,
    const IndexType nonZeroIndexes[],
    const OtherValueType nonZeroValues[],
    const Scalar zeroValue )
{
    setSameValue( n, zeroValue );
    hmemo::HArrayRef<IndexType> aNonZeroIndexes( nnz, nonZeroIndexes );
    hmemo::HArrayRef<OtherValueType> aNonZeroValues( nnz, nonZeroValues );
    fillSparseData( aNonZeroIndexes, aNonZeroValues, common::BinaryOp::COPY );
} 
  
} /* end namespace lama */

} /* end namespace scai */
