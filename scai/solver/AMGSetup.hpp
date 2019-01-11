/**
 * @file AMGSetup.hpp
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

class _AMGSetup;

// typedef std::shared_ptr<class AMGSetup> AMGSetupPtr;

typedef std::pair<common::ScalarType, std::string> AMGSetupCreateKeyType;

class COMMON_DLL_IMPORTEXPORT _AMGSetup :

    public common::Printable,
    public common::Factory<AMGSetupCreateKeyType, _AMGSetup*>
{
public:

    /**
     *  Provide a more convenient interface to the create method of the factory.
     */

    static _AMGSetup* getAMGSetup( const common::ScalarType scalarType, const std::string& setupType );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief The class AMGSetup should describe the Interace to an AMG Setup.
 *
 * @todo The current Interface of AMGSetup is just for evaluation so this should be changed to meet all requirements.
 *       (e.g. Pre and Post Smoothing)
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT AMGSetup :

    public _AMGSetup

{
public:

    AMGSetup();

    virtual ~AMGSetup();

    /**
     * @brief Create a new AMGSetup of a certain type.
     *
     */
    static AMGSetup* getAMGSetup( const std::string& setupType );

    /**
     *  Provide a more convenient interface to query for a setup type if value type is known. 
     */
    static inline bool canCreate( const std::string& setupType )
    {
        return _AMGSetup::canCreate( AMGSetupCreateKeyType( common::TypeTraits<ValueType>::stype, setupType ) );
    }

    /**
     *  Get all setup types available for this value type by using _AMGSetup::createValues 
     */
    static void getCreateValues( std::vector<std::string>& values );

    virtual void initialize( const lama::Matrix<ValueType>& coefficients ) = 0;

    virtual Solver<ValueType>& getCoarseLevelSolver() = 0;

    virtual IndexType getNumLevels() = 0;

    virtual Solver<ValueType>& getSmoother( const IndexType level ) = 0;

    virtual const lama::Matrix<ValueType>& getGalerkin( const IndexType level ) = 0;

    virtual const lama::Matrix<ValueType>& getRestriction( const IndexType level ) = 0;

    virtual const lama::Matrix<ValueType>& getInterpolation( const IndexType level ) = 0;

    virtual lama::Vector<ValueType>& getSolutionVector( const IndexType level ) = 0;

    virtual lama::Vector<ValueType>& getRhsVector( const IndexType level ) = 0;

    virtual lama::Vector<ValueType>& getTmpResVector( const IndexType level ) = 0;

    virtual std::string getCouplingPredicateInfo() const = 0;

    virtual std::string getColoringInfo() const = 0;

    virtual std::string getInterpolationInfo() const = 0;

    virtual std::string getSmootherInfo() const = 0;

    virtual std::string getCoarseLevelSolverInfo() const = 0;

    virtual void setMaxLevels( const IndexType level ) = 0;

    virtual void setMinVarsCoarseLevel( const IndexType vars ) = 0;

    virtual void setHostOnlyLevel( IndexType hostOnlyLevel );

    virtual void setHostOnlyVars( IndexType hostOnlyVars );

    virtual void setReplicatedLevel( IndexType replicatedLevel );

    virtual void setCoarseLevelSolver( SolverPtr<ValueType> solver ) = 0;

    /**
     * @brief Sets smoother for all level
     */
    virtual void setSmoother( SolverPtr<ValueType> solver ) = 0;

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
