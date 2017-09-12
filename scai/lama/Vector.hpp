/**
 * @file Vector.hpp
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
 * @brief Definition of an abstract class for distributed vectors.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distributed.hpp>

// local library
#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/Scalar.hpp>
#include <scai/lama/io/FileIO.hpp>

// others
#include <scai/common/BinaryOp.hpp>
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

class Matrix;

/** Pointer class for a vector, always use of a shared pointer. */

typedef common::shared_ptr<class Vector> VectorPtr;

/** Help class as forward declaration of enum types belonging to class Vector. */

struct _Vector
{
    /**
     * @brief VectorKind describes if a vector is dense or sparse.
     */
    typedef enum
    {
        DENSE,      //!< vector format for a dense vector
        SPARSE,     //!< vector format for a sparse vector
        JOINED,     //!< vector format for a joined vector
        UNDEFINED   //!< for convenience, always the last entry, stands also for number of entries
    } VectorKind;

    static COMMON_DLL_IMPORTEXPORT const char* kind2Str( const VectorKind vectorKind );

    static COMMON_DLL_IMPORTEXPORT VectorKind str2Kind( const char* str );

};  // struct _Vector

/** @brief Output operator<< for VectorKind prints meaningful names instead of int values */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const _Vector::VectorKind& kind );

/** Type definition for the key type used for the Vector factory.
 *
 *  The key for vector create is a pair of vector format and the value type.
 */

typedef std::pair<_Vector::VectorKind, common::scalar::ScalarType> VectorCreateKeyType;

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
class COMMON_DLL_IMPORTEXPORT Vector:

    public common::Factory<VectorCreateKeyType, Vector*>,
    public dmemo::Distributed,
    public _Vector

{
public:

    /** Help class to observe the further use of operator[] for Vector */

    class VectorElemProxy
    {
    public:

        /** Proxy constructed by ref to the array and the index value. */

        inline VectorElemProxy( Vector& vector, const IndexType i );

        /** Proxy for a vector element can be used to get its value, type conversion to Scalar
         *
         *  @returns current value of the vector element as a Scalar
         */
        inline operator Scalar() const;

        /** indexed value proxy can be assigned a scalar */

        inline VectorElemProxy& operator= ( Scalar val );

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        inline VectorElemProxy& operator= ( const VectorElemProxy& other );

    private:

        Vector& mVector;
        IndexType mIndex;
    };

    /** @brief More convenient use of the create routine of factory that avoids use of CreateKeyType.
     */
    static Vector* getVector( const VectorKind format, const common::scalar::ScalarType valueType );

    /** @brief More convenient routine to create a dense vector with certain properties.
     *
     *  @param[in] valueType specifies the type of the Vector to be created
     *  @param[in] distribution becomes the distribution of the new vector
     *  @param[in] context optional, becomes the context of the new vector
     */
    static Vector* getDenseVector(
        const common::scalar::ScalarType valueType,
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
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Vector.
     */
    typedef const Vector& ExpressionMemberType;

    /**
     * @brief Releases all allocated resources.
     */
    virtual ~Vector();

    /** Each derived vector must give info about its kind (DENSE or SPARSE). */

    virtual VectorKind getVectorKind() const = 0;

    /**
     * @brief The assignment operator assigns the result of the passed expression
     *        to this.
     *
     * The assignment operator assigns the result of the passed expression to
     * this, if necessary new memory will be allocated. The Vector will hold
     * the result of the Matrix Vector multiplication represented by
     * expression.
     *
     * @param[in] expression the input expression.
     * @return               a reference to this.
     * @throws               Exceptions thrown by the Allocator
     */

    /** this = alpha * A * x */

    Vector& operator=( const Expression_SMV& expression );

    /** this = alpha * x * A */

    Vector& operator=( const Expression_SVM& expression );

    /** this = alpha * x + beta * y */

    Vector& operator=( const Expression_SV_SV& expression );

    /** this = alpha * A * x + beta * y */

    Vector& operator=( const Expression_SMV_SV& expression );

    /** this = alpha * x * A + beta * y */

    Vector& operator=( const Expression_SVM_SV& expression );

    /** this = alpha * x */

    Vector& operator=( const Expression_SV& expression );

    /** this = alpha * x + beta */

    Vector& operator=( const Expression_SV_S& );

    /** this = x * y */

    Vector& operator=( const Expression_VV& );

    /** this = alpha * x * y */

    Vector& operator=( const Expression_SVV& );

    /** this +=  alpha * A * x */

    Vector& operator+=( const Expression_SMV& expression );

    /** this +=  alpha * x * A */

    Vector& operator+=( const Expression_SVM& expression );

    /** this +=  alpha * x */

    Vector& operator+=( const Expression_SV& expression );

    /** this -=  alpha * A * x */

    Vector& operator-=( const Expression_SMV& expression );

    Vector& operator-=( const Expression_SVM& expression );

    /** this -=  alpha * x */

    Vector& operator-=( const Expression_SV& expression );

    /**
     * @brief Assigns the values of other to the elements of this.
     *
     * @param[in] other   the vector to get values from.
     * @return            a reference to this.
     */
    Vector& operator=( const Vector& other );

    /**
     * @brief Multiplies the passed value with all elements of this.
     *
     * @param[in] value   the value to multiply all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator*=( const Scalar value );

    /**
     * @brief Multiplies the passed value with all elements of this.
     *
     * @param[in] other   the vector to multiply to do the multiplication per element
     * @return            a reference to this.
     */
    Vector& operator*=( const Vector& other );

    /**
     * @brief Divides the passed value with all elements of this.
     *
     * @param[in] value   the value to divide all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator/=( const Scalar value );

    /**
     * @brief Divide this vector by another vector element-wise
     *
     * @param[in] other   the vector to multiply to do the multiplication per element
     * @return            a reference to this.
     */
    Vector& operator/=( const Vector& other );

    /**
     * @brief Returns the addition of this and other.
     *
     * @param[in] other the vector to do the addition with.
     * @return          a reference to this.
     */
    Vector& operator+=( const Vector& other );

    /**
     * @brief Returns the subtraction of this and other.
     *
     * @param[in] other the vector to do the subtraction with.
     * @return          a reference to this.
     */
    Vector& operator-=( const Vector& other );

    /**
     * @brief Add a scalar value to all elements of this vector.
     *
     * @param[in] value   the value to add all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator+=( const Scalar value );

    /**
     * @brief Sub a scalar value to all elements of this vector.
     *
     * @param[in] value   the value to add all elements of this with.
     * @return            a reference to this.
     */
    Vector& operator-=( const Scalar value );

    /**
     * @brief Assigns the passed value to all elements of this.
     *
     * @param[in] value   the value to assign to all elements of this.
     * @return            a reference to this.
     */
    Vector& operator=( const Scalar value );

    /**
     * @brief Returns a copy of the value at the passed global index.
     *
     * @param[in] i    the global index to get the value at.
     * @return         a copy of the value at the passed global position.
     *
     * As this operator requires communication ins SPMD mode it can be very inefficient in some situations.
     */
    Scalar operator()( const IndexType i ) const
    {
        return getValue( i );
    }

    VectorElemProxy operator[]( const IndexType i )
    {
        return VectorElemProxy( *this, i );
    }

    Scalar operator[]( const IndexType i ) const
    {
        return getValue( i );
    }

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

    /**
     * @brief Sets the local values of a vector by a sparse pattern, i.e. non-zero indexes and values
     *
     * @param[in] nonZeroIndexes   array with all local indexes that have a non-zero entry
     * @param[in] nonZeroValues    array with the values for the nonZeroIndexes
     * @param[in] zeroValue        value for all indexes that do not appear in nonZeroIndexes
     *
     * The size of the both input arrays must be equal.
     *
     * Note: Implicit type conversion for the values is supported.
     */
    virtual void setSparseValues( 
        const hmemo::HArray<IndexType>& nonZeroIndexes, 
        const hmemo::_HArray& nonZeroValues,
        const Scalar zeroValue = Scalar( 0 ) ) = 0;

    /**
     * @brief Sets the local size of the vector to zero. 
     *
     * This routine has the same semantic as setValues with an empty array of any type. It is
     * a private routine as it allows a temporary inconsistency between the local part and 
     * the distribution.
     */

    virtual void clearValues() = 0;

    /**
     * This method initializes a (distributed) vector with random numbers.
     *
     * @param[in] distribution specifies the distribution of the vector
     * @param[in] fillRate for the number of non-zeros
     */
    virtual void setRandom( dmemo::DistributionPtr distribution, const float fillRate = 1.0 ) = 0;

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
     * @param[in] fileType format of the output file ("frv" for SAMG, "mtx" for MatrixMarket), default is to decide by suffix
     * @param[in] dataType representation type for output values, if set it overrides IO settings
     * @param[in] fileMode can be BINARY or FORMATTED, DEFAULT_MODE keeps default/environment settings
     */
    void writeToFile(
        const std::string& fileName,
        const std::string& fileType = "",
        const common::scalar::ScalarType dataType = common::scalar::UNKNOWN,
        const FileIO::FileMode fileMode = FileIO::DEFAULT_MODE  ) const;

    /**
     * @brief Queries the value type of the vector elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::scalar::ScalarType getValueType() const = 0;

    /**
     * @brief Returns the value at the passed global index.
     *
     * @param[in] globalIndex   the global index to get the value at.
     * @return                  a copy of the value at the passed global position.
     *
     * As this operation requires communication in SPMD mode it can be very inefficient in some situations.
     * Therefore it is recommended to query values on the local vector data with local indexes.
     */
    virtual Scalar getValue( IndexType globalIndex ) const = 0;

    /**
     *
     * @brief This methods sets/updates a value of a vector.
     *
     * Be careful: this method might throw an exception on a sparse vector, if the element is not available
     */
    virtual void setValue( const IndexType globalIndex, const Scalar value ) = 0;

    /**
     * @brief Returns the global minimum value of this.
     *
     * @return   the global minimum value of this vector.
     */
    virtual Scalar min() const = 0;

    /**
     * @brief Returns the global maximum value of this.
     *
     * @return the global maximum value of this vector.
     */
    virtual Scalar max() const = 0;

    virtual Scalar sum() const = 0;

    /**
     * @brief Returns the L1 norm of this.
     *
     * @return the L1 norm of this.
     *
     * l1Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l1Norm() const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual Scalar l2Norm() const = 0;

    /**
     * @brief Returns the max norm of this.
     *
     * @return the max norm of this.
     *
     * maxNorm computes the value of this with the largest magnitude.
     */
    virtual Scalar maxNorm() const = 0;

    /**
     * @brief Returns the max norm of the difference with another vector
     *
     *  v1.maxDiffNorm( v2 ) is equivalent to:
     *
     *  \code
     *      Vector tmp = v1 - v2;
     *      maxNorm( tmp )
     *  \endcode
     *
     *  But it avoids the temporary vector wherever possible
     */
    virtual Scalar maxDiffNorm( const Vector& other ) const = 0;

    /**
     *  Method to create a new vector of the same kind and same type
     *
     *  /code
     *    const Vector& old = ...
     *    Vector* new = Vector::create( old.getCreateValue() );
     *  /endcode
     *
     *  This routine is very important to write code that can deal with arbitrary types
     *  but does not have a template param for the value type.
     */

    virtual VectorCreateKeyType getCreateValue() const = 0;

    /**
     *  @brief Creates a new Vector of the same kind and value type, and same context
     *
     *  /code
     *    const Vector& old = ...
     *    ....
     *    Vector* new = old.newVector();
     *
     *    // is same as
     *
     *    Vector* new = Vector::create( old.getCreateValue() );
     *    new->setContextPtr( new.getContextPtr() );
     *  /endcode
     *
     *  The new vector is a zero vector, not allocated, not initialized.
     */
    virtual Vector* newVector() const = 0;

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
    virtual Vector* copy() const = 0;

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
    virtual void swap( Vector& other ) = 0;

    /** Override default implementation of Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /**
     *  @brief Assigns an arbitrary vector to this vector.
     *
     *  Common implementation for all vectors using virtual methods.
     */
    void assign( const Vector& other );

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
        const common::binary::BinaryOp op = common::binary::COPY,
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
        const common::binary::BinaryOp op = common::binary::COPY,
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
    virtual void vectorPlusVector( const Scalar& alphaS, const Vector& x, const Scalar& betaS, const Vector& y ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX * vectorY
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorTimesVector( const Scalar& alphaS, const Vector& x, const Vector& y ) = 0;

    /**
     * @brief Assignment of a 'full' vector expression vectorResult = scalarAlpha * vectorX * scalarBeta
     *
     * Each vector class has to implement its own version of this assignment. 
     */
    virtual void vectorPlusScalar( const Scalar& alphaS, const Vector& x, const Scalar& betaS ) = 0;

    /**
     * @brief Returns the dot product of this and other.
     *
     * @param[in] other   the vector to calculate the dot product with.
     * @return            the dot product of this and other
     */
    virtual Scalar dotProduct( const Vector& other ) const = 0;

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
    virtual void setVector( const Vector& other, common::binary::BinaryOp op, const bool swapArgs = false ) = 0;

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
     *      v.invert()      v.setScalar( Scalar( 1 ), binary::DIVIDE, true );
     *      v += s;         v.setScalar( s, binary::ADD, false );
     *      v *= s;         v.setScalar( s, binary::MULT, false );
     *  \endcode
     */
    virtual void setScalar( const Scalar value, common::binary::BinaryOp op, const bool swapScalar = false ) = 0;

    /**
     *  @brief Apply a unary operation for each element of the vector.
     */
    virtual void applyUnary( common::unary::UnaryOp op ) = 0;

    /**
     *  @brief Boolean reduction returns true if all elements fullfill the compare operation with a scalar.
     */
    virtual bool all( common::binary::CompareOp op, const Scalar value ) const = 0;

    /**
     *  @brief Boolean reduction returns true if elementwise comparison with other vector is true for all elements
     */
    virtual bool all( common::binary::CompareOp op, const Vector& other ) const = 0;

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
    void powBase( const Vector& other );

    /**
     *  Calculates the pow function for the vector elements with the elements of another vector.
     */
    void powExp( const Vector& other );

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
     *  Constructor of Vector for derived classes by size and/or context
     */
    explicit Vector( const IndexType size = 0, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context       is optional, will be Host context.
     */
    explicit Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Creates a copy of the passed Vector.
     *
     * @param[in] other   the Vector to take a copy from.
     *
     * Inherits size/distribution, context and content of the passed vector.
     */
    Vector( const Vector& other );

    /**
     *  @brief Swaps member variables of Vector class.
     */
    void swapVector( Vector& other );

    hmemo::ContextPtr mContext; //!< decides about location of vector operations

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** write only the local data to a file, no communication here */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const = 0;

    /** write the whole vector into a single file, can imply redistribution */

    void writeToSingleFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** same as writeLocalToFile but also communication for error handling */

    void writeToPartitionedFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
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

Vector::VectorElemProxy::VectorElemProxy( Vector& vector, const IndexType i ) :

    mVector( vector ),
    mIndex( i )

{
}

Vector::VectorElemProxy::operator Scalar() const
{
    return mVector.getValue( mIndex );
}

Vector::VectorElemProxy& Vector::VectorElemProxy::operator= ( Scalar val )
{
    mVector.setValue( mIndex, val );
    return *this;
}

Vector::VectorElemProxy& Vector::VectorElemProxy::operator= ( const Vector::VectorElemProxy& other )
{
    Scalar tmp = other.mVector.getValue( other.mIndex );
    mVector.setValue( mIndex, tmp );
    return *this;
}

IndexType Vector::size() const
{
    return getDistribution().getGlobalSize();
}

hmemo::ContextPtr Vector::getContextPtr() const
{
    return mContext;
}

} /* end namespace lama */

} /* end namespace scai */
