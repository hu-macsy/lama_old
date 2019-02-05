/**
 * @file JoinedVector.hpp
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
 * @brief Vector class that stands for the concatenation of two arbitrary vectors
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/Vector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/dmemo/JoinedDistribution.hpp>

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

    virtual void resize( scai::dmemo::DistributionPtr )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void redistribute( const scai::dmemo::RedistributePlan& )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void writeLocalToFile( FileIO& ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    virtual void readFromFile( FileIO& )
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
