/**
 * @file Criterion.hpp
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
 * @brief Criterion.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

// std
#include <ostream>
#include <memory>

namespace scai
{

namespace solver
{

template<typename ValueType> class IterativeSolver;

template<typename ValueType> class Criterion;

template<typename ValueType>
using CriterionPtr = std::shared_ptr<Criterion<ValueType> >;

/**
 * @brief Defines the operators which can be used to connect Criterions with.
 */
enum class BooleanOp
{
    AND, //!< stands for logical and composition
    OR //!< stands for logical or composition
};

/** @brief Provide output operator<< to print enum values. */

std::ostream& operator<<( std::ostream& stream, const BooleanOp op );

/**
 * @brief The class Criterion is the base class for all stopping criterions
 *        that we can use on a Solver.
 *
 * Criterion can be composed with the logical operators &&, ||, and !.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Criterion: public common::Printable
{
public:

    /**
     * @brief Creates a Criterion which evaluates to true.
     */
    Criterion();

    /**
     * @brief Creates a Criterion which evaluates to the passed boolean.
     *
     * @param[in] boolean   the boolean to set the result of this boolean condition to.
     *
     * Note: implicit type conversions of bool argument are disabled
     */
    explicit Criterion( const bool boolean );

    /**
     * @brief Creates a copy of the passed boolean condition.
     *
     * @param[in] other the boolean condition to copy.
     */
    Criterion( const Criterion& other );

    /**
     * @brief Creates a copy of the passed Criterion, which is altered by the
     *        passed modifier.
     *
     * Creates a copy of the passed Criterion, which is altered by the
     * passed modifier. The new Criterion behaves like
     * \code
     *     other.isSatisfied() != modifier
     * \endcode
     * which is equivalent to  other.isSatisfied() XOR modifier.
     *
     * @param[in] other     the boolean condition to take a altered copy from
     * @param[in] modifier  the bool value to modify the new Criterion with.
     */
    Criterion( const CriterionPtr<ValueType> other, const bool modifier );

    /**
     * @brief Creates a new Criterion which the logical composition of
     *        the passed Criterions.
     *
     * @param[in] leftChild     the left operand of the composition
     * @param[in] rightChild    the right operand of the composition
     * @param[in] operation     the logical operation to connect leftChild and
     *                          rightChild with.
     */
    Criterion( const CriterionPtr<ValueType> leftChild, const CriterionPtr<ValueType> rightChild, BooleanOp operation );

    /** Destructor. */

    virtual ~Criterion();

    /**
     * @brief Checks if this Criterion is true for the state of the passed
     *        solver.
     *
     * @param[in] solver    the solver to check the state of.
     * @return              if this Criterion is true for the state of solver.
     */
    virtual bool isSatisfied( const IterativeSolver<ValueType>& solver );

    const CriterionPtr<ValueType> getLeftChild() const;

    const CriterionPtr<ValueType> getRightChild() const;

    BooleanOp getOperation() const;

    void setLeftChild( const CriterionPtr<ValueType> leftChild );

    void setRightChild( const CriterionPtr<ValueType> rightChild );

    void setOperation( const BooleanOp operation );

    /**
     * @brief Returns true if the boolean condition has a left child
     * Do not use getLeftChild() to compare pointers. Logical operations are overwritten
     *
     * return Returns true if the boolean condition has a left child. If pointer is not set: false
     */
    bool hasLeftChild() const;

    /**
     * @brief Returns true if the boolean condition has a right child
     * Do not use getRightChild() to compare pointers. Logical operations are overwritten
     *
     * return Returns true if the boolean condition has a right child. If pointer is not set: false
     */
    bool hasRightChild() const;

    // unused: operator= always used on CriterionPtr<ValueType>
    // Criterion& operator =( const Criterion& other );

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :
    CriterionPtr<ValueType> mLeftChild;
    CriterionPtr<ValueType> mRightChild;
    BooleanOp mOperation;
    const bool mModifier; // true: negate, false: do nothing
};

template<typename ValueType>
inline CriterionPtr<ValueType> operator||( CriterionPtr<ValueType> a, CriterionPtr<ValueType> b )
{
    return CriterionPtr<ValueType>( new Criterion<ValueType>( a, b, BooleanOp::OR ) );
}

template<typename ValueType>
inline CriterionPtr<ValueType> operator&&( CriterionPtr<ValueType> a, CriterionPtr<ValueType> b )
{
    return CriterionPtr<ValueType>( new Criterion<ValueType>( a, b, BooleanOp::AND ) );
}

} /* end namespace solver */

} /* end namespace scai */
