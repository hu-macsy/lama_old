/**
 * @file DenseVector.hpp
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
 * @brief Definition of template class that stands for a dense vector of a certain type.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/_DenseVector.hpp>

// internal scai libraries
#include <scai/utilskernel/LArray.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Halo.hpp>
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/macros/throw.hpp>

// std

namespace scai
{

namespace dmemo
{
class Redistributor;
}

namespace lama
{

/**
 * @brief The template DenseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseVector:

    public _DenseVector,

    public Vector::Register<DenseVector<ValueType> >    // register at factory

{
public:

    /** Default constructor, creates empty (not initilized) vector, that is replicated (without distribution) */

    DenseVector();

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] context  the context to use for the new vector.
     */
    explicit DenseVector( hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized replicated DenseVector of the passed global size.
     *
     * @param[in] size  is the global size of the vector
     */
    explicit DenseVector( const IndexType size );

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     */
    explicit DenseVector( dmemo::DistributionPtr distribution );

    /**
     * @brief create a not initialized replicated DenseVector with a given context
     *
     * @param[in] size  is the global size of the vector
     * @param[in] context  the context to use for the new vector.
     */
    DenseVector ( const IndexType size, hmemo::ContextPtr context );

    /**
     * @brief creates a not initialized distributed DenseVector of the passed global size.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] context  the context to use for the new vector.
     */
    DenseVector ( dmemo::DistributionPtr distribution, hmemo::ContextPtr context );

    /**
     * @brief creates a replicated DenseVector of the passed size initialized to the passed value.
     *
     * @param[in] size  the size of the new DenseVector.
     * @param[in] value the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    DenseVector( const IndexType size, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief creates a distributed DenseVector of the passed global size initialized to the passed value.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] value         the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    DenseVector( dmemo::DistributionPtr distribution, const ValueType value, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * Override the default copy constructor to guarantee a deep copy.
     *
     * @param[in] other the dense vector that will be copied
     *
     * The new constructed vector has the same distribution as the input vector.
     */
    DenseVector( const DenseVector<ValueType>& other );

    /**
     * More general constructor that creates a deep copy of an arbitrary vector.
     *
     * The explicit specifier avoids implict conversions as the following example shows.
     *
     * \code
     *     subroutine sub( const DenseVector<float>& v );
     *     ...
     *     DenseVector<float> vf;
     *     DenseVector<double> vd;
     *     sub( vf );               // that is okay
     *     sub( vd );               // compile error to avoid implicit conversions
     * \endcode
     */
    explicit DenseVector( const Vector& other );

    /**
     * @brief creates a redistributed copy of the passed vector
     *
     * @param[in] other         the vector to take a copy from
     * @param[in] distribution  the distribution to use for the new vector.
     *
     * Must be valid: other.size() == distribution.getGlobalSize()
     */
    DenseVector( const Vector& other, dmemo::DistributionPtr distribution );

    /**
     * @brief creates a distributed DenseVector with given local values.
     *
     * @param[in] distribution  the distribution for the vector
     * @param[in] localValues   the local values to initialize the new DenseVector with.
     *
     * Note: localValues.size() == distribution->getLocalSize() must be valid.
     */
    DenseVector( dmemo::DistributionPtr distribution, const hmemo::_HArray& localValues );

    /**
     * @brief Constructor of a replicated DenseVector with local values
     *
     * Note: DenseVector( localValues ) is same as DenseVector( localValues, repDistribution )
     *       with repDistributon = DistributionPtr( new NoDistribution( localValues.size() )
     *
     * Usually, a replicated vector has on multiple processes the same values on each processor.
     * It might be used in a local mode where each processor has individual values but then it
     * should never be used in operations that involve global communication.
     */
    explicit DenseVector( const hmemo::_HArray& localValues );

    /**
     * @brief This constructor creates a vector with the size and values stored
     *        in the file with the given filename.
     *
     * @param[in] filename  the name of the file to read the Vector from.
     *
     * The distribution of the vector is CYCLIC(n) where n is the size of
     * the vector. So only the first processor will hold all values.
     *
     * Note: Only the first processor will read the matrix file.
     */
    explicit DenseVector( const std::string& filename );

    /**
     * @brief creates a DenseVector with the Expression alpha * x.
     *
     * @param[in] expression    alpha * x
     */
    explicit DenseVector( const Expression_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha + x.
     *
     * @param[in] expression    alpha * x + beta
     */
    explicit DenseVector( const Expression_SV_S& expression );

    /**
     *  @brief creates a DenseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    x * y
     */

    explicit DenseVector( const Expression_VV& expression );


    /**
     *  @brief creates a DenseVector with the Expression alpha * x * Y.
     *
     * @param[in] expression    alpha * x * y
     */

    explicit DenseVector( const Expression_SVV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x + beta * y.
     *
     * @param[in] expression  is alpha * x + beta * y
     */
    explicit DenseVector( const Expression_SV_SV& expression );

    /* --------------------------------------------------------------------- */

    /**
     * @brief creates a DenseVector with the Expression alpha * A * x + beta * y.
     *
     * @param[in] expression     alpha * A * x + beta * y
     */
    explicit DenseVector( const Expression_SMV_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x * A + beta * y.
     *
     * @param[in] expression     alpha * x * A + beta * y
     */
    explicit DenseVector( const Expression_SVM_SV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * A * x.
     *
     * @param[in] expression     alpha * A * x
     */
    explicit DenseVector( const Expression_SMV& expression );

    /**
     * @brief creates a DenseVector with the Expression alpha * x * A.
     *
     * @param[in] expression     alpha * x * A
     */
    explicit DenseVector( const Expression_SVM& expression );

    /**
     * @brief creates a DenseVector with the Expression A * x.
     *
     * @param[in] expression     A * x
     */
    explicit DenseVector( const Expression_MV& expression );

    /**
     * @brief creates a DenseVector with the Expression x * A.
     *
     * @param[in] expression     x * A
     */
    explicit DenseVector( const Expression_VM& expression );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~DenseVector();

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /** Implememenation of pure routine Vector::allocate. */

    virtual void allocate( const IndexType n );

    /** Allocate and initialize this vector with raw values 
     *
     *  @tparam OtherValueType data type of the raw data
     *  @param[in] size becomes the size of the vector and specifies number of entries in values
     *  @param[in] values is the array with the raw data
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

    /** Override the default assignment operator.  */

    DenseVector& operator=( const DenseVector<ValueType>& other );

    /** Reimplement Vector::operator= as otherwise a constructor of DenseVector might be called */

    DenseVector& operator=( const Scalar );

    // All other assignment operators are inherited from class Vector, but using is required

    using Vector::operator=;
    using Vector::assign;

    /**
     * Implementation of pure method Vector::fillRandom 
     */
    virtual void fillRandom( const IndexType bound );

    /** Implementation of pure method Vector::fillSparseRandom */

    virtual void fillSparseRandom( const float fillRate, const IndexType bound );

    /** Implementation of pure method _DenseVector::fillRange */

    virtual void fillRange( const Scalar startValue, const Scalar inc );

    /** Sort all elements of this vector.
     *
     *  Currently, sorting is only possible on block distributed vectors.
     *  Keep in mind that this operation might introduce a new (general block) distribution
     *  for the vector.
     *
     *  @param[out] perm         permutation vector with the global indexes of original positions
     *  @param[in]  ascending    flag if sorting is ascending or descending
     *
     *  Note:  oldVector[ perm ] ->  newVector
     */

    virtual void sort( DenseVector<IndexType>& perm, bool ascending );

    /** Same sort but without perm vector. */

    virtual void sort( bool ascending );

    /** Checking whether all values in a dense vector are sorted.
     *
     *  This method can only be applied for block distributions.
     */
    virtual bool isSorted( bool ascending ) const;

    /** Compute a global prefix sum for the elements of this vector. 
     *
     *  \code
     *   A = { 1, 2, 3, 4, 5, 6 }
     *   A.scan();
     *   A = { 1, 3, 6, 10, 15, 21 }
     *  \endcode  
     *
     *  Note: This vector must be block distributed, otherwise it throws an exception
     */
    void scan();

    /** Same as scan but it uses another input vector. */

    void scan( const DenseVector<ValueType>& other );

    /** Implementation of Vector::getValueType */

    virtual common::scalar::ScalarType getValueType() const;

    /**
     * Implementation of pure method Vector::setDenseValues.
     */
    virtual void setDenseValues( const hmemo::_HArray& values );

    /**
     * Implementation of pure method Vector::fillSparseData
     */
    virtual void fillSparseData( 
        const hmemo::HArray<IndexType>& nonZeroIndexes,
        const hmemo::_HArray& nonZeroValues,
        const common::binary::BinaryOp op );

    /**
     * Implementation of Vector::copy with covariant return type.
     */
    virtual DenseVector* copy() const;

    /**
     * Implementation of Vector::newVector with covariant return type.
     */
    virtual DenseVector* newVector() const;

    /**
     * @brief get a non constant reference to local values of this Dense Vector.
     *
     * @return  a non constant reference to the local values of this.
     */

    inline utilskernel::LArray<ValueType>& getLocalValues();

    /**
     * @brief get a constant reference to local values of this Dense Vector.
     *
     * @return  a constant reference to the local values of this.
     */
    inline const utilskernel::LArray<ValueType>& getLocalValues() const;

    /**
     * @brief Get a reference to the halo temp array of this Dense Vector.
     *
     * @return  a reference to the halo values of this.
     *
     * The halo array can be used as temporary array to keep values of neighbored
     * processors. It avoids reallocation of memory for the values.
     */

    inline utilskernel::LArray<ValueType>& getHaloValues() const;

    virtual Scalar getValue( IndexType globalIndex ) const;

    virtual void setValue( const IndexType globalIndex, const Scalar value );

    virtual Scalar min() const;

    virtual Scalar max() const;

    virtual Scalar sum() const;

    virtual Scalar l1Norm() const;

    virtual Scalar l2Norm() const;

    virtual Scalar maxNorm() const;

    /** Implementation of pure method Vector::maxDiffNorm */

    virtual Scalar maxDiffNorm( const Vector& other ) const;

    /** Implementation of pure method Vector::all */

    virtual bool all( common::binary::CompareOp op, const Scalar value ) const;

    /** Implementation of pure method Vector::all */

    virtual bool all( common::binary::CompareOp op, const Vector& other ) const;

    virtual void swap( Vector& other );

    /** Reset a dense vector with a new array of local values and a new distribution.
     *
     *  @param[in,out] newValues array with the new values before, contains old values afterwards
     *  @param[in]     newDist is the new distribution
     *
     *  This methods throws an exception if the size of the new local values does not fit with the
     *  local size of the new distribution.
     *
     *  This method is more efficient than calling a constructor for DenseVector if the array
     *  newValues is no more needed afterwards.
     */
    void swap( hmemo::HArray<ValueType>& newValues, dmemo::DistributionPtr newDist );

    virtual void writeAt( std::ostream& stream ) const;

    /** Implementation of pure method Vector::vectorPlusVector */

    virtual void vectorPlusVector( const Scalar& alphaS, const Vector& x, const Scalar& betaS, const Vector& y );

    /** vectorPlusVector with one zero term */

    void assignScaledVector( const Scalar& alpha, const Vector& x );

    /** vectorPlusVector with aliased */

    void axpy( const Scalar& alpha, const Vector& x );

    /** Implementation of pure method Vector::vectorTimesVector */

    virtual void vectorTimesVector( const Scalar& alphaS, const Vector& x, const Vector& y );

    /** Implementation of pure method Vector::vectorPlusScalar */

    virtual void vectorPlusScalar( const Scalar& alphaS, const Vector& x, const Scalar& betaS );

    /** Assign this vector with a scalar values, does not change size, distribution. */

    virtual void assign( const Scalar value );

    /** Implementation of pure method Vector::setScalar */

    virtual void setScalar( const Scalar value, common::binary::BinaryOp op, const bool swapScalar = false );

    /** Implementation of pure method Vector::setVector */

    virtual void setVector( const Vector& other, const common::binary::BinaryOp op, const bool swapScalar = false );

    /** Implementation of pure method Vector::applyUnary */

    virtual void applyUnary( common::unary::UnaryOp op );

    /** Setting this vector by gathering vector elements from another vector.
     *
     *  @param[in] source is the vector from which elements are gathered
     *  @param[in] index  are the elements needed from other processors
     *  @param[in] op     specifies how to combine elements with existing ones
     *
     *  If op is COPY, this vector will have the same size and distribution as index, otherwise
     *  this vector and the other vector must have the same size and distribution.
     */
    virtual void gather(
        const DenseVector<ValueType>& source,
        const DenseVector<IndexType>& index,
        const common::binary::BinaryOp op = common::binary::COPY );

    /** Scattering values from another vector into this vector
     *
     *  @param[in] index  specifies positions where to update values
     *  @param[in] source values that are scattered
     *  @param[in] op     specifies how to combine elements with existing ones
     *
     *  *this[ index ] = source, index and source must have same size/distribution
     */
    virtual void scatter(
        const DenseVector<IndexType>& index,
        const DenseVector<ValueType>& source,
        const common::binary::BinaryOp op = common::binary::COPY );

    /**
     * @brief Implementation of pure method Vector::buildLocalValues.
     */
    virtual void buildLocalValues( 
        hmemo::_HArray& localValues,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /**
     * @brief Implementation of pure method Vector::gatherLocalValues.
     *
     * For a dense vector this is just a corresponding gather of the HArray containing
     * all local values.
     */
    virtual void gatherLocalValues(
        hmemo::_HArray& localValues,
        const hmemo::HArray<IndexType>& localIndexes,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    virtual Scalar dotProduct( const Vector& other ) const;

    using Vector::prefetch; // prefetch() with no arguments

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( dmemo::DistributionPtr distribution );

    /** Implementation of pure method Vector::redistribute */

    virtual void redistribute( const dmemo::Redistributor& redistributor );

protected:

    using Vector::mContext;

private:

    /** Allocate local values array with correct size at right context */

    void allocate();

    /** Static method for sorting of DenseVector */

    static void sortImpl(
        DenseVector<IndexType>* perm,
        DenseVector<ValueType>* out,
        DenseVector<ValueType>& in,
        bool descending );

    utilskernel::LArray<ValueType> mLocalValues; //!< my local values of vector

    /** array that might be used to keep halo values of vector, avoids reallocation of memory for halo values */

    mutable utilskernel::LArray<ValueType> mHaloValues;

    /** Implementation of Vector::writeLocalToFile */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::scalar::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** Implementation of Vector::readLocalFromFile */

    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = nIndex );

    /** Implementation of Vector::clearValues */

    virtual void clearValues();

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();

    virtual VectorCreateKeyType getCreateValue() const;
};

/* ------------------------------------------------------------------------- */

template<typename ValueType>
utilskernel::LArray<ValueType>& DenseVector<ValueType>::getLocalValues()
{
    return mLocalValues;
}

template<typename ValueType>
const utilskernel::LArray<ValueType>& DenseVector<ValueType>::getLocalValues() const
{
    return mLocalValues;
}

template<typename ValueType>
utilskernel::LArray<ValueType>& DenseVector<ValueType>::getHaloValues() const
{
    return mHaloValues;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::setRawData( const IndexType size, const OtherValueType values[] )
{
    allocate( size );
    // use LAMA array reference to avoid copy of the raw data
    hmemo::HArrayRef<OtherValueType> valuesArrayRef( size, values );
    // use mContext instead of context to avoid NULL pointer
    utilskernel::HArrayUtils::assign( mLocalValues, valuesArrayRef, mContext );
}

} /* end namespace lama */

} /* end namespace scai */
