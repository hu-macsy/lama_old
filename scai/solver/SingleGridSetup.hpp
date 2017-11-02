/**
 * @file SingleGridSetup.hpp
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

class SingleGridSetup:

    public AMGSetup,
    public AMGSetup::Register<SingleGridSetup>    // register at factory

{
public:
    SingleGridSetup();
    virtual ~SingleGridSetup();

    virtual void initialize( const lama::_Matrix& coefficients );

    virtual Solver& getCoarseLevelSolver();

    virtual unsigned int getNumLevels();

    virtual Solver& getSmoother( const unsigned int level );

    virtual const lama::_Matrix& getGalerkin( const unsigned int level );

    virtual const lama::_Matrix& getRestriction( const unsigned int level );

    virtual const lama::_Matrix& getInterpolation( const unsigned int level );

    virtual lama::_Vector& getSolutionVector( const unsigned int level );

    virtual lama::_Vector& getRhsVector( const unsigned int level );

    virtual lama::_Vector& getTmpResVector( const unsigned int level );

    virtual std::string getCouplingPredicateInfo() const;

    virtual std::string getColoringInfo() const;

    virtual std::string getInterpolationInfo() const;

    virtual std::string getSmootherInfo() const;

    virtual std::string getCoarseLevelSolverInfo() const;

    virtual void setCoarseLevelSolver( SolverPtr solver );

    /**
     * @brief Usually sets smoother for all level, for this case overwrites the coarse level solver
     */
    virtual void setSmoother( SolverPtr solver );

    // just a dummy function
    virtual void setMaxLevels( const unsigned int )
    {
    }
    // just a dummy function
    virtual void setMinVarsCoarseLevel( const unsigned int )
    {
    }

    // Get the key used for registration in factory

    static std::string createValue();

    // Create routine used for factory

    static AMGSetup* create();

private:

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    SolverPtr    mSolver;

    std::unique_ptr<lama::_Matrix> mIdentity;

    std::unique_ptr<lama::_Vector> mSolutionVector;
    std::unique_ptr<lama::_Vector> mRhsVector;
    std::unique_ptr<lama::_Vector> mTmpResVector;

};

} /* end namespace solver */

} /* end namespace scai */
