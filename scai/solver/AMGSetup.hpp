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
 * @author Thomas Brandes, Jiri Kraus
 * @date 28.10.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/solver/Solver.hpp>

// internal scai libraries
#include <scai/common/Factory.hpp>

namespace scai
{

namespace solver
{

/** 
 *  @brief Untyped base class for all typed AMG classes
 *
 *  This common base class provides a common factory and a common logger.
 */
class _AMGSetup;

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

    /**
     *  Initalization of the AMG setup by the 'square' input matrix.
     *
     *  @param[in] mainSystemMatrix is the input matrix used for the AMG setup
     *
     *  Note: row and column distribution of the matrix should be equal.
     *
     *  The context and format of the input matrix is inherited to all galerkin
     *  and interpolation/restriction matrices.
     */
    void initialize( const lama::Matrix<ValueType>& mainSystemMatrix );

    /**
     *  @brief Get the number of levels
     *
     *  Note: this method returns 0 if not initialized, 1 if there is only the main matrix used.
     */
    IndexType getNumLevels() const;

    /**
     *  @brief Getter for the galerkin matrix on a certain level
     * 
     *  @param[in] level must be valid index, $0 \le level \lt getNumLevels()$
     * 
     *  The matrix on level 0 is the matrix that has been used for initialization. The
     *  higher the level, the less entries the matrix has.
     */
    const lama::Matrix<ValueType>& getGalerkin( const IndexType level );

    /**
     *  @brief Getter for the restriction matrix to restrict form one level to the 
     */
    const lama::Matrix<ValueType>& getRestriction( const IndexType level );

    /**
     *  @brief Getter for the interpolation matrix to interpolate a matrix at a level from the next higher level.
     */
    const lama::Matrix<ValueType>& getInterpolation( const IndexType level );

    lama::Vector<ValueType>& getSolutionVector( const IndexType level );

    lama::Vector<ValueType>& getRhsVector( const IndexType level );

    lama::Vector<ValueType>& getTmpResVector( const IndexType level );

    /** 
     *  @brief Creation of interpolation, galerkin and interpolation matrices must
     *         be provided by each derived class.
     *
     *  Note: Derived classes should add the matrices for the next level by calling 
     *        the protected method AMGSetup::addNextLevel.
     */ 
    virtual void createMatrixHierarchy() = 0;

    /**
     *  @brief Get a solver to be used as a smoother.
     *
     *  This method must be implemented by all derived classes.
     */
    virtual SolverPtr<ValueType> createSolver( bool isCoarseLevel ) = 0;

    virtual std::string getCouplingPredicateInfo() const = 0;

    virtual std::string getColoringInfo() const = 0;

    virtual std::string getInterpolationInfo() const = 0;

    /**
     *  @brief Set the maximal number of levels
     *
     *  By default, there is no restriction for the maximal number of levels.
     */
    void setMaxLevels( const IndexType level );

    /**
     *  @brief Set the minimal size for a matrix that it should have on the coarsest level
     *
     *  This size is used as stopping criterion for creation of a next level.
     */
    void setMinVarsCoarseLevel( const IndexType vars );

    void setHostOnlyLevel( IndexType hostOnlyLevel );

    void setHostOnlyVars( IndexType hostOnlyVars );

    void setReplicatedLevel( IndexType replicatedLevel );

    /**
     * @brief Sets default smoother that is used on each level
     *
     * Note: this method should be called before initialize.
     */
    void setSmoother( SolverPtr<ValueType> solver );

    /**
     * @brief Sets default solver to be used on the coarsest level
     *
     * Note: this method should be called before initialize.
     */
    void setCoarseLevelSolver( SolverPtr<ValueType> solver );

    /**
     * @brief Getter routine for the solver on a level.
     */
    Solver<ValueType>& getSmoother( const IndexType level );

    /**
     * @brief Getter routine for the solver on the coarsest level.
     */
    Solver<ValueType>& getCoarseLevelSolver();

    std::string getSmootherInfo() const;

    std::string getCoarseLevelSolverInfo() const;

protected:

    IndexType mHostOnlyLevel;   //<! determines how many of the coarsest grids are kept on Host

    IndexType mHostOnlyVars;    //<! matrices with a size less or equal have always host context

    IndexType mReplicatedLevel;  //<! determines how many of the coarses levels will be replicated

    IndexType mMinVarsCoarseLevel;  //!< no further level if matrix has less equal variables

    IndexType mMaxLevels;           //!< maximal number of matrices in the hiearchy 

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    /** 
     *  This method shoud be called by the derived AMG setup class to set the matrices for the next level.
     *
     *  @param[in] interpolationMatrix is the 
     *
     *  This call takes over the ownership of all matrices. The matrices will not be changed, i.e.
     *  the calling method might use e.g. a reference of the galerkin matrix to compute the next level.
     */
    void addNextLevel( 
        std::unique_ptr<lama::Matrix<ValueType>> interpolationMatrix,
        std::unique_ptr<lama::Matrix<ValueType>> galerkinMatrix = std::unique_ptr<lama::Matrix<ValueType>>(),
        std::unique_ptr<lama::Matrix<ValueType>> restrictionMatrix = std::unique_ptr<lama::Matrix<ValueType>>() );

private:

    /**
     *  Help routine that converts the matrices (galerkin, restriction, interpolation) on the different levels
     */
    void convertMatrixHierarchy();

    /**
     *  Help routine that sets up the vectors used on each level.
     */
    void createVectorHierarchy();

    /**
     *  Help routine that sets up the solvers/smoothers used on each level.
     */
    void createSolverHierarchy();

    /** @brief pointer reference to the main coefficient matrix
     *
     *  It is assumed that the matrix is not modified during the whole lifetime of the setup.
     */
    const lama::Matrix<ValueType>* mMainGalerkinMatrix;   

    /** 
     *  @brief Data structure that keeps the galerkin matrices on all levels, except the main one.
     *
     *  The main system matrix is not managed here as we do not have the ownership.
     */
    std::vector<std::unique_ptr<lama::Matrix<ValueType>>> mGalerkinMatrices;

    /** 
     *  @brief Data structure that keeps the interpolation matrices between all levels 
     */
    std::vector<std::unique_ptr<lama::Matrix<ValueType>>> mInterpolationMatrices;

    /** 
     *  @brief Data structure that keeps the restriction matrices between all levels 
     *
     *  Note: the restriction matrix is the tranpose of the interpolation matrix.
     */
    std::vector<std::unique_ptr<lama::Matrix<ValueType>>> mRestrictionMatrices;

    /**
     *  @brief Setup keeps the allocated data for the solution, rhs, residual vectors
     */
    std::vector<lama::DenseVector<ValueType>> mSolutionHierarchy;
    std::vector<lama::DenseVector<ValueType>> mRhsHierarchy;
    std::vector<lama::DenseVector<ValueType>> mTmpResHierarchy;

    /**
     *  @brief Data structure that keeps the solver on each level 
     */
    std::vector<SolverPtr<ValueType> > mSolverHierarchy;

    /**
     *  @brief Default solver that should be used on the coarsest level.
     */
    SolverPtr<ValueType> mCoarseLevelSolver;

    /**
     *  @brief Default solver to be used for smoothing on each level.
     */
    SolverPtr<ValueType> mSmoother; 
};

} /* end namespace solver */

} /* end namespace scai */
