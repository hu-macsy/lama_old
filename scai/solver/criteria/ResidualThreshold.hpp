/**
 * @file ResidualThreshold.hpp
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
 * @brief ResidualThreshold.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/criteria/Criterion.hpp>

#include <scai/lama/norm/Norm.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The possible check modes for the residue threshold
 */
enum class ResidualCheck
{
    Absolute, //!< stands for an absolute residue reduction
    Relative, //!< stands for an relative residue reduction
    Divergence //!< stands for an relative residue increase
};

/** @brief Provide output operator<< to print enum values. */

std::ostream& operator<<( std::ostream& stream, const ResidualCheck checkMode );

/**
 * @brief ResidualThreshold is a stopping criterion of a solver which checks the
 *        residue of a solver for certain criteria.
 *
 * ResidualThreshold is a stopping criterion of a solver which checks the residue
 * of a solver to fall below a configured threshold. The threshold is either checked
 * absolute or relative to the initial residue.
 *
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ResidualThreshold: public Criterion<ValueType>
{
public:
    /*
     * @brief Empty constructor first. Configured by setter methods.
     */
    ResidualThreshold();

    /**
     * @brief Creates a RedisualThreshold stopping criteria from the passed norm.
     *
     * @param[in] norm      the norm to use for residue calculation
     */
    ResidualThreshold( const lama::NormPtr<ValueType> norm );

    /**
     * @brief Creates a RedisualThreshold stopping criteria from the passed parameters.
     *
     * @param[in] norm      the norm to use for residue calculation
     * @param[in] precision the desired precision
     * @param[in] checkMode if the residue should be check for an absolute or a relative
     *                      reduction.
     */
    ResidualThreshold( const lama::NormPtr<ValueType> norm, const ValueType precision, const ResidualCheck checkMode );

    /**
     * @brief Creates a copy of the passed ResidualThreshold.
     *
     * @param[in] other the ResidualThreshold to take a copy from.
     */
    ResidualThreshold( const ResidualThreshold& other );

    /**
     * @brief Destroys this ResidualThreshold.
     */
    virtual ~ResidualThreshold();

    /** Implementation of pure method Criterion<ValueType>::copy with covariant return type. */

    virtual ResidualThreshold<ValueType>* copy() const;

    bool isSatisfied( const IterativeSolver<ValueType>& solver );

    ResidualCheck getCheckMode() const;
    IndexType getIterationExtrema() const;
    ValueType getFirstNormResult() const;
    const lama::NormPtr<ValueType> getNorm() const;
    ValueType getPrecision() const;

    void setCheckMode( ResidualCheck mCheckMode );
    void setFirstNormResult( ValueType firstNormResult );
    void setPrecision( ValueType precision );

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_USING( Criterion<ValueType>::logger );

private:

    const lama::NormPtr<ValueType> mNorm;

    ResidualCheck mCheckMode;
    RealType<ValueType> mPrecision;
    RealType<ValueType> mFirstNormResult;
};

} /* end namespace solver */

} /* end namespace scai */
