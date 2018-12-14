/**
 * @file DenseVector.hpp
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
#include <scai/dmemo/Distribution.hpp>
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


// forward declaration required as we do not include SparseVector.hpp here

template<typename ValueType>
class SparseVector;

/**
 * @brief The template class DenseVector represents a distributed 1D Vector with elements of type ValueType.
 *
 * @tparam ValueType the value type for the vector values.
 *
 * In contrarary to a sparse vector a dense vector has allocated memory for each entry. If the vector
 * is distributed on multiple processors, each processor has allocated only memory for its local part, 
 * i.e. for the elements that are owned by it.
 *
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
    using _Vector::getDistributionPtr;
    using _Vector::size;
    using _Vector::assign;
    using _Vector::prefetch; // prefetch() with no arguments


    using Vector<ValueType>::operator=;
    using Vector<ValueType>::getValueType;
    using Vector<ValueType>::binaryOp;    

    /** Default constructor, creates empty (not initilized) vector, that is replicated (without distribution) */

    explicit DenseVector( hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief Create a replicated DenseVector of a certain size initialized with same value for all elements.
     *
     * @param[in] size  the size of the new DenseVector.
     * @param[in] value the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     *
     * Attention: DEPRECATED.
     *
     *  \code
     *   DenseVector<ValueType> v( 100, 1 );
     *   auto v = fill<DenseVector<ValueType>>( 100, 1 );
     *  \endcode
     */
    DenseVector( const IndexType size, const ValueType value, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief creates a distributed DenseVector of the passed global size initialized to the passed value.
     *
     * @param[in] distribution  the distribution to use for the new vector.
     * @param[in] value         the value to assign to all elements of the new DenseVector.
     * @param[in] context   specifies optionally the context where dense vector should reside
     */
    DenseVector( dmemo::DistributionPtr distribution, const ValueType value, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * Override the default copy constructor to guarantee a deep copy.
     *
     * @param[in] other the dense vector that will be copied
     *
     * The new constructed vector has the same distribution as the input vector.
     */
    DenseVector( const DenseVector<ValueType>& other );

    /**
     * @brief Override the default move constructor with appropriate version.
     *
     * @param[in] other the dense vector from which data is moved
     *
     * The attribute noexcept is mandatory to make the use of dense vectors possible in C++ container classes.
     */
    DenseVector( DenseVector<ValueType>&& other ) noexcept;

    /**
     * @brief creates a distributed DenseVector with given local values.
     *
     * @param[in] distribution  the distribution for the vector
     * @param[in] localValues   the local values to initialize the new DenseVector with.
     * @param[in] context       becomes the context for memory/operations on this vector.
     *
     * Note: localValues.size() == distribution->getLocalSize() must be valid.
     */
    DenseVector( 
        dmemo::DistributionPtr distribution, 
        hmemo::HArray<ValueType> localValues, 
        hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief Constructor of a replicated DenseVector with global values
     *
     *  \code
     *     HArray<double> values( { 1, 2, 3, 4, 5 } );
     *     DenseVector<double>( values );
     *     // short form of: DenseVector<double>( std::make_shared<NoDistriubiton<values.size()>, values );
     *  \endcode
     *
     * Usually, a replicated vector has the same values on each processor.
     * It might be used in a local mode where each processor has individual values but then you
     * should exactly know what you do.
     */
    explicit DenseVector( hmemo::HArray<ValueType> localValues, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief releases all allocated resources.
     */
    virtual ~DenseVector();

    /* ==================================================================== */
    /*   split up member variables                                          */
    /* ==================================================================== */

    void splitUp( dmemo::DistributionPtr& dist, hmemo::HArray<ValueType>& localvalues );

    /** 
     * @brief Implementation of pure method _Vector::getVectorKind() 
     */

    inline virtual VectorKind getVectorKind() const;

    /**
     * @brief Implementation of pure method _Vector::isConsistent 
     */
    virtual bool isConsistent() const;

    /* --------------------------------------------------------------------- */

    /**
     *  @brief Implememenation of pure routine _Vector::allocate. 
     *
     *  Please keep in mind that this method does not initialize the values 
     *  of the dense vector. So it should be used only in combination with
     *  methods that set the values afterwards.
     *
     *  \code
     *     DenseVector<double> v;       DenseVector<double> v( 100, 1.0 );     DenseVector<double> v;
     *     v.allocate( 100 );                                                  v.setSameValue( 100, 1.0 );
     *     v = 1.0;             
     *  \endcode
     */
    virtual void allocate( const IndexType n );

    /** Implememenation of pure routine _Vector::allocate. */

    virtual void allocate( dmemo::DistributionPtr distribution );

    /* --------------------------------------------------------------------- */

    /** Override the default assignment operator.  */

    DenseVector& operator=( const DenseVector<ValueType>& other );

    /**
     * @brief Move assignment operator for dense vector.
     *
     * @param[in] other the dense vector to move
     * @return reference to this array for further assignments
     *
     * The other dense vector will be a zero vector afterwards.
     */
    DenseVector<ValueType>& operator=( DenseVector<ValueType>&& other ) noexcept;

    // All other assignment operators are inherited from class Vector, but using is required

    /** Implementation of pure method _Vector::asign 
     *
     *  uses metaprogramming to call assignImpl with actual type and kind
     */
    virtual void assign( const _Vector& other );

    /** 
     *  @brief Specific implementation of assign
     */
    template<typename OtherValueType>
    void assignImpl( const Vector<OtherValueType>& other );

    /** 
     *  @brief Specific implementation of assign for a sparse vector of a certain type.
     */
    template<typename OtherValueType>
    void assignSparse( const SparseVector<OtherValueType>& other );

    /** 
     *  @brief Specific implementation of assign for a dense vector of a certain type.
     */
    template<typename OtherValueType>
    void assignDense( const DenseVector<OtherValueType>& other );

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

    inline hmemo::HArray<ValueType>& getLocalValues();

    /**
     * @brief get a constant reference to local values of this Dense Vector.
     *
     * @return  a constant reference to the local values of this.
     */
    inline const hmemo::HArray<ValueType>& getLocalValues() const;

    /**
     * @brief Get a reference to the halo temp array of this Dense Vector.
     *
     * @return  a reference to the halo values of this.
     *
     * The halo array can be used as temporary array to keep values of neighbored
     * processors. It avoids reallocation of memory for the values.
     */

    inline hmemo::HArray<ValueType>& getHaloValues() const;

    /** @brief Implementation of pure method Vector<ValueType>::getValue */

    virtual ValueType getValue( IndexType globalIndex ) const;

    /** @brief Implementation of pure method Vector<ValueType>::setValue */

    virtual void setValue( const IndexType globalIndex, const ValueType value );

    virtual ValueType min() const;

    virtual ValueType max() const;

    virtual ValueType sum() const;

    /** Implementation of pure method Vector<ValueType>::l1Norm() for dense vector */

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation of pure method Vector<ValueType>::l2Norm() for dense vector */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation of pure method Vector<ValueType>::maxNorm() for dense vector */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation of pure method Vector<ValueType>::maxDiffNorm */

    virtual RealType<ValueType> maxDiffNorm( const Vector<ValueType>& other ) const;

    /** Implementation of pure method Vector<ValueType>::all */

    virtual bool all( common::CompareOp op, const ValueType value ) const;

    /** Implementation of pure method Vector<ValueType>::all */

    virtual bool all( common::CompareOp op, const Vector<ValueType>& other ) const;

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

    /** Implementation of pure method Vector<ValueType>::unaryOp for dense vector. */

    void unaryOp( const Vector<ValueType>& x, common::UnaryOp op );

    /** Implementation of pure method Vector<ValueType>::binaryOp for dense vector. */

    void binaryOp( const Vector<ValueType>& x, common::BinaryOp op, const Vector<ValueType>& y );

    /** Implementation of pure method Vector<ValueType>::selectComplexPart for dense vector. */

    virtual void selectComplexPart( Vector<RealType<ValueType> >& x, const common::ComplexPart kind ) const;

    /** Implementation of pure method Vector<ValueType>::buildComplex for dense vector. */

    virtual void buildComplex( const Vector<RealType<ValueType> >& x, const Vector<RealType<ValueType> >& y );

    /** Implementation of pure method Vector<ValueType>::binaryOpScalar for dense vector. */

    void binaryOpScalar( const Vector<ValueType>& x, const ValueType& alpha, const common::BinaryOp op, const bool swap );

    /** Implementation of pure method Vector<ValueType>::vectorPlusVector */

    virtual void vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta, const Vector<ValueType>& y );

    /** vectorPlusVector with aliased */

    void axpy( const ValueType& alpha, const Vector<ValueType>& x );

    /** Implementation of pure method Vector<ValueType>::vectorTimesVector */

    virtual void vectorTimesVector( const ValueType& alpha, const Vector<ValueType>& x, const Vector<ValueType>& y );

    /** Implementation of pure method Vector<ValueType>::vectorPlusScalar */

    virtual void vectorPlusScalar( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta );

    /** Implementation for pure method Vector<ValueType>::setScalar */

    virtual void setScalar( const ValueType& alpha );

    /** Implementation of pure method _Vector::setVector */

    virtual void setVector( const _Vector& other, const common::BinaryOp op, const bool swapScalar = false );

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

    /** Implementation of pure method Vector<ValueType>::dotProduct */

    virtual ValueType dotProduct( const Vector<ValueType>& other ) const;

    virtual void prefetch( const hmemo::ContextPtr location ) const;

    virtual void wait() const;

    virtual size_t getMemoryUsage() const;

    virtual void redistribute( dmemo::DistributionPtr distribution );

    /** Implementation of pure method _Vector::redistribute */

    virtual void redistribute( const dmemo::Redistributor& redistributor );

    /** Implementation of pure method _Vector::resize */

    virtual void resize( const dmemo::DistributionPtr distribution );

private:

    using _Vector::setDistributionPtr;

    /** Allocate local values array with correct size at right context */

    void allocate();

    /** Static method for sorting of DenseVector */

    static void sortImpl(
        DenseVector<IndexType>* perm,
        DenseVector<ValueType>* out,
        DenseVector<ValueType>& in,
        bool descending );

    hmemo::HArray<ValueType> mLocalValues; //!< my local values of vector

    /** array that might be used to keep halo values of vector, avoids reallocation of memory for halo values */

    mutable hmemo::HArray<ValueType> mHaloValues;

    /** Implementation of _Vector::writeLocalToFile */

    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileIO::FileMode fileMode ) const;

    /** Implementation of _Vector::readLocalFromFile */

    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = invalidIndex );

    /** Implementation of _Vector::clearValues */

    virtual void clearValues();

public:

    // static methods, variables to register create routine in Vector factory of base class.

    static _Vector* create();

    // key for factory

    static VectorCreateKeyType createValue();
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
DenseVector<ValueType> linearDenseVector( 
    const IndexType n,
    const ValueType startValue, 
    const ValueType inc, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
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
DenseVector<ValueType> linearDenseVector( 
    dmemo::DistributionPtr distribution, 
    const ValueType startValue, 
    const ValueType inc, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
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
hmemo::HArray<ValueType>& DenseVector<ValueType>::getLocalValues()
{
    return mLocalValues;
}

template<typename ValueType>
const hmemo::HArray<ValueType>& DenseVector<ValueType>::getLocalValues() const
{
    return mLocalValues;
}

template<typename ValueType>
hmemo::HArray<ValueType>& DenseVector<ValueType>::getHaloValues() const
{
    return mHaloValues;
}

/* ------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
