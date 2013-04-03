/**
 * @file SolverFactory.hpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief SolverFactory.hpp
 * @author kbuschulte
 * @date 02.06.2012
 * $Id$
 */
#ifndef LAMA_SOLVERFACTORY_HPP_
#define LAMA_SOLVERFACTORY_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/solver/creator/SolverCreator.hpp>
#include <lama/solver/Solver.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

namespace qi = boost::spirit::qi;

/**
 * registry for Solver
 * provides interface for getting a new CG solver
 */
class LAMA_DLL_IMPORTEXPORT SolverFactory
{
public:
    typedef qi::symbols<char,SolverCreator::RuleType> TypeToCreatorMap;
    typedef qi::symbols<char,SolverPtr> SolverInstanceMap;

    /**
     * @brief Returns the static factory object to access the factory
     *
     * @return  static factory object
     */
    static SolverFactory& getFactory();

    virtual ~SolverFactory();

    /**
     * @brief Add/register a SolverCreator to the Factory
     *
     * @param[in] type             The SolverType which is used to create a solver instance
     * @param[in] SolverCreator    The creator to create solver instances
     */
    void addSolverCreator( const std::string& type, SolverCreator::RuleType& SolverCreator );

    const TypeToCreatorMap& getCreatorRuleSymbols();

    /**
     * @brief Checks if a solver instance is registered
     *
     * @param[in] solverName    The name/id of the registered solver
     *
     * @return Returns true if the solver is available/registered
     */
    bool hasSolver( const std::string& solverName ) const;

    /**
     * @brief Returns a registered solver instance
     *
     * @param[in] solverName    The name/id of the registered solver
     *
     * @return Returns the solver if the solver is available/registered
     */
    SolverPtr getSolver( const std::string& solverName );

    /**
     * @brief Returns the solver instance map
     *
     * @return Returns the solver instance map
     */
    const SolverInstanceMap& getSolverInstanceMap();

    /**
     * @brief Add/register a solver instance with a given name/id to the Factory
     *
     * @param[in] solver    The solver instance which will be added to the registry/map
     */
    void addSolver( Solver* solver );

    /**
     *  Releases all registered solvers and solver creators.
     */
    static void release();

private:
    /**
     * @brief Default constructor of Solver Factory
     * Create and manage a collection of solvers
     */
    SolverFactory();

    /* solver types holding all possible solver managers to create a solver */
    TypeToCreatorMap mCreatorRuleMap;

    /* solver types holding all possible solver managers to create a solver */
    SolverInstanceMap mSolverInstanceMap;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

#endif // LAMA_SOLVERFACTORY_HPP_
