/**
 * @file SingleGridSetup.hpp
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
 * @brief SingleGridSetup.hpp
 * @author Jiri Kraus
 * @date 27.10.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/AMGSetup.hpp>

// local library
#include <scai/solver/Solver.hpp>

namespace scai
{

namespace solver
{

/** 
 *  @brief Trivial implementation for an AMG setup.
 *
 *  This AMG setup has only 1 level that contains the system matrix itself.
 *  As default, it takes a Jacobi solver (10 iterations steps) for smoothing
 */
template<typename ValueType>
class SingleGridSetup:

    public AMGSetup<ValueType>,
    public _AMGSetup::Register<SingleGridSetup<ValueType> >   // register at factory

{
public:

    SingleGridSetup();

    virtual ~SingleGridSetup();

    virtual std::string getCouplingPredicateInfo() const;

    virtual std::string getColoringInfo() const;

    virtual std::string getInterpolationInfo() const;

    // Get the key used for registration in factory

    static AMGSetupCreateKeyType createValue();

    // Create routine used for factory

    static _AMGSetup* create();

private:

    /**
     *   Implementation of pure method AMGSetup::createMatrixHieararchy()
     */
    virtual void createMatrixHierarchy();

    /**
     *   @brief Implementation of pure method AMGSetup::createSolver
     */
    virtual SolverPtr<ValueType> createSolver( bool isCoarseLevel );

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
