/**
 * @file Criterion.cpp
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
 * @brief Criterion.cpp
 * @author Kai Buschulte
 * @date 25.07.2011
 */

// hpp
#include <scai/solver/criteria/Criterion.hpp>

// local library
#include <scai/solver/IterativeSolver.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

std::ostream& operator<<( std::ostream& stream, const BooleanOp op )
{
    switch ( op )
    {
        case BooleanOp::AND :
            stream << "AND";
            break;
        case BooleanOp::OR :
            stream << "OR";
            break;
        default:
            stream << "<unknown_boolean_op>";
    }

    return stream;
}

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, Criterion<ValueType>::logger, "Criterion" );

template<typename ValueType>
Criterion<ValueType>::Criterion() : 

    mOperation( BooleanOp::AND ), 
    mModifier( true )
{
}

template<typename ValueType>
Criterion<ValueType>::Criterion( const bool boolean ) : 

    mOperation( BooleanOp::AND ), 
    mModifier( boolean )
{
}

template<typename ValueType>
Criterion<ValueType>::Criterion( const Criterion<ValueType>& other ) : 

    Printable( other ), 
    mLeftChild( other.mLeftChild ), 
    mRightChild( other.mRightChild ), 
    mOperation( other.mOperation ), 
    mModifier( false )
{
}

template<typename ValueType>
Criterion<ValueType>::Criterion( const CriterionPtr<ValueType> other, const bool modifier ) : 

    mLeftChild( other->mLeftChild ), 
    mRightChild( other->mRightChild ), 
    mOperation( other->mOperation ), 
    mModifier( modifier == other->mModifier )
{
}

template<typename ValueType>
Criterion<ValueType>::Criterion( const CriterionPtr<ValueType> leftChild, const CriterionPtr<ValueType> rightChild, BooleanOp operation ) : 

    mLeftChild( leftChild ), 
    mRightChild( rightChild ), 
    mOperation( operation ), 
    mModifier( false )
{
    SCAI_ASSERT_DEBUG( leftChild, "Left child is a NULL pointer." );
    SCAI_ASSERT_DEBUG( rightChild, "Right child is a NULL pointer." );
}

template<typename ValueType>
Criterion<ValueType>::~Criterion()
{
}

template<typename ValueType>
bool Criterion<ValueType>::isSatisfied( const IterativeSolver<ValueType>& solver )
{
    SCAI_LOG_INFO( logger, "isSatisfied: " << *this );
    bool satisfied = true;

    if ( !hasLeftChild() && !hasRightChild() )
    {
        return mModifier;
    }

    if ( !hasLeftChild() )
    {
        satisfied = mRightChild->isSatisfied( solver );
    }
    else if ( !hasRightChild() )
    {
        satisfied = mLeftChild->isSatisfied( solver );
    }
    else if ( mOperation == BooleanOp::AND )
    {
        satisfied = mLeftChild->isSatisfied( solver ) && mRightChild->isSatisfied( solver );
    }
    else
    {
        satisfied = mLeftChild->isSatisfied( solver ) || mRightChild->isSatisfied( solver );
    }

    return mModifier != satisfied; //XOR
}

template<typename ValueType>
const CriterionPtr<ValueType> Criterion<ValueType>::getLeftChild() const
{
    return mLeftChild;
}

template<typename ValueType>
const CriterionPtr<ValueType> Criterion<ValueType>::getRightChild() const
{
    return mRightChild;
}

template<typename ValueType>
BooleanOp Criterion<ValueType>::getOperation() const
{
    return mOperation;
}

template<typename ValueType>
void Criterion<ValueType>::setLeftChild( const CriterionPtr<ValueType> leftChild )
{
    mLeftChild = leftChild;
}

template<typename ValueType>
void Criterion<ValueType>::setRightChild( const CriterionPtr<ValueType> rightChild )
{
    mRightChild = rightChild;
}

template<typename ValueType>
bool Criterion<ValueType>::hasLeftChild() const
{
    return mLeftChild.get();
}

template<typename ValueType>
bool Criterion<ValueType>::hasRightChild() const
{
    return mRightChild.get();
}

template<typename ValueType>
void Criterion<ValueType>::setOperation( const BooleanOp operation )
{
    mOperation = operation;
}

template<typename ValueType>
void Criterion<ValueType>::writeAt( std::ostream& stream ) const
{
    if ( hasLeftChild() && hasRightChild() ) //boolean operation
    {
        if ( !mModifier )
        {
            stream << "!";
        }

        stream << "(" << *getLeftChild();

        if ( getOperation() == BooleanOp::AND )
        {
            stream << " && ";
        }
        else
        {
            stream << " || ";
        }

        stream << *getRightChild() << ")";
    }
    else if ( hasLeftChild() || hasRightChild() )
    {
        if ( !mModifier )
        {
            stream << "!";
        }

        stream << "(";

        if ( hasLeftChild() )
        {
            stream << *getLeftChild();
        }
        else
        {
            stream << *getRightChild();
        }

        stream << ")";
    }
    else //leaf
    {
        stream << "leaf<" << mModifier << ">";
    }
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Criterion, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
