/**
 * @file JoinedVector.hpp
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
 * @brief Vector class that stands for the concatenation of two arbitrary vectors
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/Vector.hpp>
#include <scai/lama/Scalar.hpp>

#include "JoinedDistribution.hpp"

namespace scai

{

namespace lama

{

/** An object of this class encapsulates a joined vector [ V1 ; V2 ]
 */

template<typename ValueType>
class JoinedVector : public Vector<ValueType>
{

public:

    JoinedVector ( Vector<ValueType>& v1, Vector<ValueType>& v2 ) :

        Vector<ValueType>( dmemo::DistributionPtr( new dmemo::JoinedDistribution( v1.getDistributionPtr(), v2.getDistributionPtr() ) ) ),

        mV1( v1 ),
        mV2( v2 )

    {
        SCAI_ASSERT_EQ_ERROR( v1.getValueType(), v2.getValueType(), "joined vectors must have same type" )
    }

    JoinedVector ( Vector<ValueType>* v1, Vector<ValueType>* v2 ) :

        Vector<ValueType>( dmemo::DistributionPtr( new dmemo::JoinedDistribution( v1->getDistributionPtr(), v2->getDistributionPtr() ) ) ),

        mV1( *v1 ),
        mV2( *v2 )

    {
        SCAI_ASSERT_EQ_ERROR( v1->getValueType(), v2->getValueType(), "joined vectors must have same type" )
        pV1.reset( v1 );
        pV2.reset( v2 );
    }

    ~JoinedVector()
    {
    }

    void assign( const _Vector& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    Vector<ValueType>& first()
    {
        return mV1;
    }

    Vector<ValueType>& second()
    {
        return mV2;
    }

    const Vector<ValueType>& first() const
    {
        return mV1;
    }

    const Vector<ValueType>& second() const
    {
        return mV2;
    }

    virtual bool isConsistent() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual VectorKind getVectorKind() const
    {
        return VectorKind::JOINED;
    }

    virtual void setDenseValues( const scai::hmemo::_HArray& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void fillSparseData( const scai::hmemo::HArray<int>&, const scai::hmemo::_HArray&, common::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void clearValues()
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void fillRandom( IndexType )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void fillSparseRandom( float, IndexType )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual ValueType getValue( IndexType ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void setValue( IndexType, ValueType )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual ValueType min() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual ValueType max() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual ValueType sum() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual RealType<ValueType> l1Norm() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual RealType<ValueType> l2Norm() const
    {
        ValueType res = dotProduct( *this );
        return common::Math::sqrt( res );
    }

    virtual RealType<ValueType> maxNorm() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual RealType<ValueType> maxDiffNorm( const Vector<ValueType>& ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual VectorCreateKeyType getCreateValue() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** This method might be used in solvers for temporary vectors. */

    virtual Vector<ValueType>* newVector() const
    {
        return new JoinedVector( mV1.newVector(), mV2.newVector() );
    }

    virtual Vector<ValueType>* copy() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void swap( _Vector& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void buildLocalValues( scai::hmemo::_HArray&, scai::common::BinaryOp, scai::hmemo::ContextPtr ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void gatherLocalValues( scai::hmemo::_HArray&, const scai::hmemo::HArray<int>&, scai::common::BinaryOp, scai::hmemo::ContextPtr ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method Vector<ValueType>::setScalar */

    virtual void setScalar( const ValueType& s )
    {
        mV1.setScalar( s );
        mV2.setScalar( s );
    }

    virtual void vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& v1, const ValueType& beta, const Vector<ValueType>& v2 )
    {
        SCAI_LOG_INFO( logger, "joinedVector = " << alpha << " * v1 + " << beta << " * v2" )

        SCAI_ASSERT_EQ_ERROR( v1.getDistribution(), this->getDistribution(), "distribution mismatch" );
        SCAI_ASSERT_EQ_ERROR( v2.getDistribution(), this->getDistribution(), "distribution mismatch" );

        const JoinedVector& jV1 = reinterpret_cast<const JoinedVector&>( v1 );
        const JoinedVector& jV2 = reinterpret_cast<const JoinedVector&>( v2 );

        mV1.vectorPlusVector( alpha, jV1.first(), beta, jV2.first() );
        mV2.vectorPlusVector( alpha, jV1.second(), beta, jV2.second() );
    }

    virtual void vectorTimesVector( const ValueType& alpha, const Vector<ValueType>& v1, const Vector<ValueType>& v2 )
    {
        SCAI_LOG_INFO( logger, "joinedVector = " << alpha << " * v1 * v2, v1 = " << v1 << ", v2 = " << v2 )

        SCAI_ASSERT_EQ_ERROR( v1.getDistribution(), this->getDistribution(), "distribution mismatch" );
        SCAI_ASSERT_EQ_ERROR( v2.getDistribution(), this->getDistribution(), "distribution mismatch" );

        const JoinedVector& jV1 = reinterpret_cast<const JoinedVector&>( v1 );
        const JoinedVector& jV2 = reinterpret_cast<const JoinedVector&>( v2 );

        mV1.vectorTimesVector( alpha, jV1.first(), jV2.first() );
        mV2.vectorTimesVector( alpha, jV1.second(), jV2.second() );
    }

    virtual void vectorPlusScalar( const ValueType&, const Vector<ValueType>&, const ValueType& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual ValueType dotProduct( const Vector<ValueType>& other ) const
    {
        SCAI_ASSERT_EQ_ERROR( other.getDistribution(), this->getDistribution(), "distribution mismatch" );

        const JoinedVector& jOther = reinterpret_cast<const JoinedVector&>( other );

        ValueType s1 = mV1.dotProduct( jOther.first() );
        ValueType s2 = mV2.dotProduct( jOther.second() );

        return s1 + s2;
    }

    virtual void setVector( const _Vector&, scai::common::BinaryOp, bool )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void unaryOp( const Vector<ValueType>&, const common::UnaryOp )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void binaryOp( const Vector<ValueType>&, const common::BinaryOp, const Vector<ValueType>& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void binaryOpScalar( const Vector<ValueType>&, const ValueType&, const common::BinaryOp, bool )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual bool all( scai::common::CompareOp, ValueType ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual bool all( scai::common::CompareOp, const Vector<ValueType>& ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void prefetch( scai::hmemo::ContextPtr ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void wait() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual size_t getMemoryUsage() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void allocate( scai::dmemo::DistributionPtr )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void allocate( IndexType )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void redistribute( scai::dmemo::DistributionPtr )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void redistribute( const scai::dmemo::Redistributor& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void writeLocalToFile( const std::string&, const std::string&, scai::common::ScalarType, FileIO::FileMode ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual IndexType readLocalFromFile( const std::string&, IndexType, IndexType )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Override default implementation of Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "JoinedVector( " << mV1 << ", " << mV2 << " )";
    }

    virtual void concatenate( dmemo::DistributionPtr, const std::vector<const Vector<ValueType> *>& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method Vector<ValueType>::selectComplexPart for dense vector. */

    virtual void selectComplexPart( Vector<RealType<ValueType> >&, const common::ComplexPart ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void buildComplex( const Vector<RealType<ValueType> >&, const Vector<RealType<ValueType> >& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

private:

    using _Vector::logger;

    Vector<ValueType>& mV1;
    Vector<ValueType>& mV2;

    // joined vector might be owner of the individual parts

    VectorPtr<ValueType> pV1;
    VectorPtr<ValueType> pV2;
};

}

}
