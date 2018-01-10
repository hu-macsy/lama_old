/**
 * @file ResidualStagnation.hpp
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
 * @brief ResidualStagnation.hpp
 * @author Kai Buschulte
 * @date 25.07.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/criteria/Criterion.hpp>

// internal scai libraries
#include <scai/lama/norm/Norm.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType> class IterativeSolver;

/**
 * @brief Stagnation stopping criterion.
 *
 * @param IterativeSolver The type solver in which this criterion is used.
 *        Most commonly this is the abstract sblas::IterativeSolver type.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ResidualStagnation: public Criterion<ValueType>
{

public:

    /**
     * @brief Creates a new Stagnation based stopping criterion.
     * TODO[doxy] param in or out?
     *
     * @param      norm The norm which shall be used to check the criterion.
     *             Has to be compatible with the _Matrix/Vector types of
     *             the solver.
     */
    ResidualStagnation( lama::NormPtr<ValueType> norm );

    /**
     * @brief Creates a new Stagnation based stopping criterion.
     * TODO[doxy] params in or out?
     *
     * @param      norm The norm which shall be used to check the criterion.
     *             Has to be compatible with the _Matrix/Vector types of
     *             the solver.
     * @param      lookback The amount of residual-norm-calculation results used
     *             for the criterion check.
     * @param      precision TODO[doxy] Complete Description.
     */
    ResidualStagnation( lama::NormPtr<ValueType> norm, IndexType lookback, ValueType precision );

    /**
     * @brief Creates a copy of the passed ResidualStagnation object.
     *
     * @param[in] other   ResidualStagnation object to be copied.
     */
    ResidualStagnation( const ResidualStagnation& other );

    /** Destructor. */

    virtual ~ResidualStagnation();

    /**
     * @brief Implement pure method but return with the covariant return type
     */
    virtual ResidualStagnation<ValueType>* copy() const;

    /**
     * @brief Checks if the criterion is satisfied.
     *
     * The Stagnation criterion needs a certain amount (specified by
     * m_lookback) to be able to perform the check. If this method is
     * called and not enough data has been gathered, false is returned.
     * Anyhow the norm of the residual is calculated and stored. This
     * process continues, until enough data has been gathered.
     *
     * @param[in] solver   The solver which calls this method. Needed for
     *                     residual computation.
     * @return   Boolean value which determines if the criterion is
     *           satisfied or not.
     */
    virtual bool isSatisfied( const IterativeSolver<ValueType>& solver );

    /**
     * @brief Gets the amount of calculations used for the criterion check.
     *
     * @return   Amount of calculations.
     */
    IndexType getLookback() const;

    /**
     * @brief Gets the Norm which the criterion uses for the check.
     *
     * @return   The Norm used.
     */
    const lama::NormPtr<ValueType> getNorm() const;

    /**
     * @brief Returns the precision for floating point number calculations related to the Norm.
     *
     * @return   The precision given in a Scalar object.
     */
    const ValueType getPrecision() const;

    /**
     * @brief Sets the amount of calculations used for the criterion check.
     *
     * @param[in] lookback   TODO[doxy] Complete Description.
     */
    void setLookback( const IndexType lookback );

    /**
     * @brief Sets the precision for floating point number calculations related to the Norm.
     *
     * @param[in] precision   should be real value
     */
    void setPrecision( const ValueType precision );

    virtual void writeAt( std::ostream& stream ) const;

private:

    /**
     * @brief The norm used in the residual-norm-calculation.
     */
    const lama::NormPtr<ValueType> mNorm;

    /**
     * @brief The amount of calculations used for the criterion check.
     */
    IndexType mLookback;

    /**
     * @brief Stores the results of the residual-norm-calculations.
     */
    std::vector<RealType<ValueType> > mLastResidualNorms;

    /**
     * @brief Index. Needed for circular use of the vector.
     */
    IndexType mNextEntry;

    /**
     * @brief Used to check, if enough data have been gathered to check the
     *        criterion.
     */
    bool mEntriesReady;

    /**
     * @brief The precision used for the stagnation check.
     */
    RealType<ValueType> mPrecision;
};

} /* end namespace solver */

} /* end namespace scai */
