/**
 * @file AMGSetup.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief AMGSetup.hpp
 * @author Jiri Kraus
 * @date 28.10.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/solver/Solver.hpp>

// internal scai libraries
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/common/Factory.hpp>

namespace scai
{

namespace solver
{

typedef common::shared_ptr<class AMGSetup> AMGSetupPtr;

/**
 * @brief The class AMGSetup should describe the Interace to an AMG Setup.
 *
 * @todo The current Interface of AMGSetup is just for evaluation so this should be changed to meet all requirements.
 *       (e.g. Pre and Post Smoothing)
 */
class COMMON_DLL_IMPORTEXPORT AMGSetup :

    public common::Factory<std::string, AMGSetup*>,
    public common::Printable

{
public:

    AMGSetup();

    virtual ~AMGSetup();

    virtual void initialize( const lama::Matrix& coefficients ) = 0;

    virtual Solver& getCoarseLevelSolver() = 0;

    virtual unsigned int getNumLevels() = 0;

    virtual Solver& getSmoother( const unsigned int level ) = 0;

    virtual const lama::Matrix& getGalerkin( const unsigned int level ) = 0;

    virtual const lama::Matrix& getRestriction( const unsigned int level ) = 0;

    virtual const lama::Matrix& getInterpolation( const unsigned int level ) = 0;

    virtual lama::Vector& getSolutionVector( const unsigned int level ) = 0;

    virtual lama::Vector& getRhsVector( const unsigned int level ) = 0;

    virtual lama::Vector& getTmpResVector( const unsigned int level ) = 0;

    virtual std::string getCouplingPredicateInfo() const = 0;

    virtual std::string getColoringInfo() const = 0;

    virtual std::string getInterpolationInfo() const = 0;

    virtual std::string getSmootherInfo() const = 0;

    virtual std::string getCoarseLevelSolverInfo() const = 0;

    virtual void setMaxLevels( const unsigned int level ) = 0;

    virtual void setMinVarsCoarseLevel( const unsigned int vars ) = 0;

    virtual void setHostOnlyLevel( IndexType hostOnlyLevel );

    virtual void setHostOnlyVars( IndexType hostOnlyVars );

    virtual void setReplicatedLevel( IndexType replicatedLevel );

    virtual void setCoarseLevelSolver( SolverPtr solver ) = 0;

    /**
     * @brief Sets smoother for all level
     */
    virtual void setSmoother( SolverPtr solver ) = 0;

protected:

    IndexType mHostOnlyLevel;

    IndexType mHostOnlyVars;

    IndexType mReplicatedLevel;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;
};

} /* end namespace solver */

} /* end namespace scai */
