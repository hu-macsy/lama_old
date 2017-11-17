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
#include <scai/lama/Vector.hpp>

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


// forward declaration required as we do not include DenseVector.hpp here

template<typename ValueType>
class SparseVector;

/**
 * @brief The template DenseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseVector:

    public Vector<ValueType>,
    public _Vector::Register<DenseVector<ValueType> >    // register at factory

{
public:

    using _Vector::logger;
    using _Vector::getDistribution;
    using _Vector::getContextPtr;
    using _Vector::getContext;
    using _Vector::readFromFile;
    using _Vector::setDistributionPtr;
    using _Vector::getDistributionPtr;
    using _Vector::size;
    using _Vector::assign;

    using Vector<ValueType>::operator=;
    using Vector<ValueType>::getValueType;

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
    explicit DenseVector( const _Vector& other );

    /**
     * @brief creates a redistributed copy of the passed vector
     *
     * @param[in] other         the vector to take a copy from
     * @param[in] distribution  the distribution to use for the new vector.
     *
     * Must be valid: other.size() == distribution.getGlobalSize()
     */
    DenseVector( const _Vector& other, dmemo::DistributionPtr distribution );

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

    /* --------------------------------------------------------------------- */

    /** 
     * @brief Implementation of pure method _Vector::getVectorKind() 
     */

    inline virtual VectorKind getVectorKind() const;

    /**
     * @brief Implementation of pure method _Vector::isConsistent 
     */
    virtual bool isConsistent() const;

    /* --------------------------------------------------------------------- */

    /** Implememenation of pure routine _Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /** Implememenation of pure routine _Vector::allocate. */

    virtual void allocate( const IndexType n );

    /** Override the default assignment operator.  */

    DenseVector& operator=( const DenseVector<ValueType>& other );

    /** Reimplement _Vector::operator= as otherwise a constructor of DenseVector might be called */

    DenseVector& operator=( const Scalar );

    // All other assignment operators are inherited from class Vector, but using is required

    /** Implementation of pure method _Vector::asign 
     *
     *  uses metaprogramming to call assignImpl with actual type and kind
     */
    virtual void assign( const _Vector& other );

    /** 
     *  @brief Specific implementation of assign for a sparse vector of a certain type.
     */
    template<typename OtherValueType>
    void assignImpl( const SparseVector<OtherValueType>& other );

    /** 
     *  @brief Specific implementation of assign for a dense vector of a certain type.
     */
    template<typename OtherValueType>
    void assignImpl( const DenseVector<OtherValueType>& other );

    /**
     * Implementation of pure method _Vector::fillRandom 
     */
    virtual void fillRandom( const IndexType bound );

    /** Implementation of pure method _Vector::fillSparseRandom */

    virtual void fillSparseRandom( const float fillRate, const IndexType bound );

    /**
     * @brief This method fills/initializes an allocated vector 
     *
     * @param[in] startValue value for the first element
     * @param[in] inc increment between the elements
     *
     */
    void fillLinearValues( const ValueType startValue, const ValueType inc );

    /** Implemenation of pure method _Vector::cat */

    virtual void concatenate( dmemo::DistributionPtr dist, const std::vector<const _Vector*>& vectors );

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

    /**
     * Implementation of pure method _Vector::setDenseValues.
     */
    virtual void setDenseValues( const hmemo::_HArray& values );

    /**
     * Implementation of pure method _Vector::fillSparseData
     */
    virtual void fillSparseData( 
        const hmemo::HArray<IndexType>& nonZeroIndexes,
        const hmemo::_HArray& nonZeroValues,
        const common::BinaryOp op );

    /**
     * Implementation of _Vector::copy with covariant return type.
     */
    virtual DenseVector* copy() const;

    /**
     * Implementation of _Vector::newVector with covariant return type.
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

    /** @brief Implementation of pure method Vector<ValueType>::getValue */

    virtual ValueType getValue( IndexType globalIndex ) const;

    /** @brief Implementation of pure method Vector<ValueType>::setValue */

    virtual void setValue( const IndexType globalIndex, const ValueType value );

    virtual ValueType min() const;

    virtual ValueType max() const;

    virtual ValueType sum() const;

    virtual NormType<ValueType> l1Norm() const;

    virtual NormType<ValueType> l2Norm() const;

    virtual NormType<ValueType> maxNorm() const;

    /** Implementation of pure method _Vector::maxDiffNorm */

    virtual NormType<ValueType> maxDiffNorm( const _Vector& other ) const;

    /** Implementation of pure method _Vector::all */

    virtual bool all( common::CompareOp op, const Scalar value ) const;

    /** Implementation of pure method _Vector::all */

    virtual bool all( common::CompareOp op, const _Vector& other ) const;

    virtual void swap( _Vector& other );

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

    /** Implementation of pure method _Vector::vectorPlusVector */

    virtual void vectorPlusVector( const Scalar& alphaS, const _Vector& x, const Scalar& betaS, const _Vector& y );

    /** vectorPlusVector with one zero term */

    void assignScaledVector( const Scalar& alpha, const _Vector& x );

    /** vectorPlusVector with aliased */

    void axpy( const Scalar& alpha, const _Vector& x );

    /** Implementation of pure method _Vector::vectorTimesVector */

    virtual void vectorTimesVector( const Scalar& alphaS, const _Vector& x, const _Vector& y );

    /** Implementation of pure method _Vector::vectorPlusScalar */

    virtual void vectorPlusScalar( const Scalar& alphaS, const _Vector& x, const Scalar& betaS );

    /** Assign this vector with a scalar values, does not change size, distribution. */

    virtual void assign( const Scalar value );

    /** Implementation of pure method _Vector::setScalar */

    virtual void setScalar( const Scalar value, common::BinaryOp op, const bool swapScalar = false );

    /** Implementation of pure method _Vector::setVector */

    virtual void setVector( const _Vector& other, const common::BinaryOp op, const bool swapScalar = false );

    /** Implementation of pure method _Vector::applyUnary */

    virtual void applyUnary( common::UnaryOp op );

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
        const common::BinaryOp op = common::BinaryOp::COPY );

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
        const common::BinaryOp op = common::BinaryOp::COPY );

    /**
     * @brief Implementation of pure method _Vector::buildLocalValues.
     */
    virtual void buildLocalValues( 
        hmemo::_HArray& localValues,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    /**
     * @brief Implementation of pure method _Vector::gatherLocalValues.
     *
     * For a dense vector this is just a corresponding gather of the HArray containing
     * all local values.
     */
    virtual void gatherLocalValues(
        hmemo::_HArray& localValues,
        const hmemo::HArray<IndexType>& localIndexes,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() ) const;

    virtual ValueType dotProduct( const _Vector& other ) const;

    using _Vector::prefetch; // prefetch() with no arguments

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( dmemo::DistributionPtr distribution );

    /** Implementation of pure method _Vector::redistribute */

    virtual void redistribute( const dmemo::Redistributor& redistributor );

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

    /** Implementation of _Vector::writeLocalToFile */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** Implementation of _Vector::readLocalFromFile */

    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = nIndex );

    /** Implementation of _Vector::clearValues */

    virtual void clearValues();

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static _Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();

    virtual VectorCreateKeyType getCreateValue() const;
};

/* ------------------------------------------------------------------------- */
/*  Functions that return a DenseVector                                      */
/* ------------------------------------------------------------------------- */

/**
 * @brief create a replicated vector and fill it with linear values
 *
 * @param[in] n becomes the size of the vector
 * @param[in] startValue value for the first element, $vector[0] = startValue$
 * @param[in] inc increment between two elements, i.e. $vector[i+1] - vector[i] = inc$
 * @param[in] ctx specifies the context of the vector.
 *
 * Note: if ctx is specified, the vector will be allocated in its memory.
 */
template<typename ValueType>
DenseVector<ValueType> linearValuesVector( 
    const IndexType n,
    const ValueType startValue, 
    const ValueType inc, 
    hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    DenseVector<ValueType> result( ctx );
    result.allocate( n );
    result.fillLinearValues( startValue, inc );
    return result;
}
/**
 * @brief create a distributed vector and fill it with linear values
 *
 * @param[in] distribution determines global/local size of the vector
 * @param[in] startValue value for the first elemen
 * @param[in] inc increment between the element
 * @param[in] ctx context where the vector is allocated
 */
template<typename ValueType>
DenseVector<ValueType> linearValuesVector( 
    dmemo::DistributionPtr distribution, 
    const ValueType startValue, 
    const ValueType inc, 
    hmemo::ContextPtr ctx = hmemo::ContextPtr() )
{
    DenseVector<ValueType> result( ctx );
    result.allocate( distribution );
    result.fillLinearValues( startValue, inc );
    return result;
}

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
VectorKind DenseVector<ValueType>::getVectorKind() const
{
    return VectorKind::DENSE;
}

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

} /* end namespace lama */

} /* end namespace scai */
