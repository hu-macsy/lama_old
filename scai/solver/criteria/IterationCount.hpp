/**
 * @file IterationCount.hpp
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
 * @brief IterationCount.hpp
 * @author Kai Buschulte
 * @date 21.07.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/criteria/Criterion.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType> class IterativeSolver;

/**
 * @brief IterationCount is a stopping criterion of a solver which checks the
 *        number of iterations a solver has executed.
 *
 * IterationCount is a stopping criterion of a solver which checks the number of
 * iterations a solver has executed. IterationCount is either true if the number
 * of iterations is larger or smaller than a configured number of iterations.
 * Depending on the configured IterationCheckMode.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT IterationCount: public Criterion<ValueType>
{
public:
    /**
     * @brief Creates a IterationCount with iterationExtream 1.
     */
    IterationCount();

    /**
     * @brief Creates a IterationCount with the passed iterationExtrema.
     *
     * @param[in] iterationExtrema   the number of iterations a solver should execute at most or minimal
     */
    IterationCount( const IndexType iterationExtrema );

    /**
     * @brief Creates a copy of the passed IterationCount object.
     *
     * @param[in] other   IterationCount object to be copied.
     */
    IterationCount( const IterationCount& other );

    /** Destructor. */

    virtual ~IterationCount();

    /**
     * @brief Apply this citerion for a given solver
     */
    virtual bool isSatisfied( const IterativeSolver<ValueType>& solver );

    /**
     * @brief Getter of the iteration extrema.
     *
     * @return   the iteration extrema.
     */
    IndexType getIterationExtrema() const;

    /**
     * @brief Setter of the iteration extrema.
     *
     * @param[in] iterationExtrema   the new iteration extrema.
     */
    void setIterationExtrema( IndexType iterationExtrema );

    virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    IndexType mIterationExtrema;
};

} /* end namespace solver */

} /* end namespace scai */
