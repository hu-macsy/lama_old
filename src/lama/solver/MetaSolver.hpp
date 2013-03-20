/**
 * @file MetaSolver.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief MetaSolver.hpp
 * @author Kai Buschulte
 * @date 07.05.2012
 * $Id$
 */

#ifndef LAMA_METASOLVER_HPP_
#define LAMA_METASOLVER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/Solver.hpp>

// others
#include <lama/solver/creator/SolverCreator.hpp>
#include <lama/solver/creator/SolverFactory.hpp>

// spirit
#include <boost/spirit/include/qi.hpp>

#include <string>
#include <fstream>

namespace lama
{

namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class LAMA_DLL_IMPORTEXPORT MetaSolver: public Solver
{
public:
    typedef std::string::const_iterator StringIterator;

    /**
     * @brief Creates a solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    MetaSolver( const std::string& id );

    /**
     * @brief Create a solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    MetaSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Create a solver with a given ID and a default-logger. The solver configuration is
     * defined by a std::string.
     *
     * @param id            The ID of the solver.
     * @param configuration Configuration for solver(s)
     */
    MetaSolver( const std::string& id, const std::string& configuration );

    /**
     * @brief Create a solver with a given ID and a given logger. The solver configuration is
     * defined by a std::string.
     *
     * @param id            The ID of the solver.
     * @param configuration Configuration for solver(s)
     * @param logger        The logger which shall be used by the solver
     */
    MetaSolver( const std::string& id, const std::string& configuration, LoggerPtr logger );

    MetaSolver( const MetaSolver& other );

    /**
     * @brief virtual destructor of MetaSolver
     */
    virtual ~MetaSolver();

    /**
     * @brief Used to initialize a solver with a certain matrix A
     *        from A*u=f.
     *
     * This method initializes a solver with a certain coefficient-matrix.
     * The only thing it does is storing the matrix pointer as a member
     * for derived solver classes to use it. The caller delegates the
     * property of the pointer to the Solver instance.
     *
     * This method may be overwritten by base classes which desire
     * more complex initialization.
     *
     * @param coefficients The matrix A from A*u=f.
     */
    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Used to initialize the preconditioner of a solver
     *        with a different matrix A'
     *
     * This method initializes a solver with a certain coefficient-matrix.
     * The only thing it does is storing the matrix pointer as a member
     * for derived solver classes to use it. The caller delegates the
     * property of the pointer to the Solver instance.
     *
     * This method may be overwritten by base classes which desire
     * more complex initialization.
     *
     * @param[in] coefficients   The matrix A from A*u=f.
     * @param[in] level          TODO[doxy] Complete Description.
     */
    virtual void initializePreconditioner( const Matrix& coefficients, LogLevel::LogLevel level );

    /**
     * @brief Solves the equation system based on the given rhs.
     *
     * The solver needs to be initialized first with the matrix from
     * the equation to solve, e.g.
     *     A from A*u=f (call solver::initialize(A) for example)
     * This method is abstract. It has to be implemented by a class which inherits
     * from this class.
     */
    virtual void solveImpl();

    /*
     * @brief If arg is a file it reads the files content. If not it calls the direct interpretation
     *
     * @param[in] arg   Input string which is a file name or a configuration string
     */
    void interpreteArgument( const std::string arg );

    /*
     * @brief Parses the given configuration and constructs the MetaSolver
     *
     * Parses the given configuration in mConfiguration and constructs the MetaSolver with the
     * described settings.
     *
     * @param[in] configuration   TODO[doxy] Complete Description.
     */
    void parseConfiguration( const std::string& configuration );

    struct MetaSolverRuntime: SolverRuntime
    {
        MetaSolverRuntime();
        virtual ~MetaSolverRuntime();

        SolverPtr mRootSolver;
    };

    virtual const MetaSolverRuntime& getConstRuntime() const;

private:
    virtual SolverPtr copy();

    virtual MetaSolverRuntime& getRuntime();

    MetaSolverRuntime mMetaSovlerRuntime;

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

};

class SolverConfigGrammar: public qi::grammar<std::string::const_iterator,void(),ascii::space_type>
{
public:
    typedef std::string::const_iterator StringIterator;
    SolverConfigGrammar();

    void setRootSolver( SolverPtr solver );

    SolverPtr getRootSolver();

private:
    typedef SolverCreator::RuleType RuleType;

    qi::rule<StringIterator,void(),qi::locals<RuleType>,ascii::space_type> mRSolverConfiguration;

    qi::rule<StringIterator,void(),ascii::space_type> mRCriteriaConfiguration;

    qi::rule<StringIterator,void(),ascii::space_type> mRLoggerConfiguration;

    qi::rule<StringIterator,void(),ascii::space_type> mRConfiguration;

    qi::rule<StringIterator,void(),ascii::space_type> mRConfigurationSequence;

    SolverPtr mRootSolver;

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

};

//template <typename Iterator>
//struct SkipGrammar : qi::grammar<Iterator, void()>
//{
//    SkipGrammar();
//
//    qi::rule<Iterator, void()> comment;
//    qi::rule<Iterator, void()> start;
//};

}//namespace lama

#endif /* METASOLVER_HPP_ */
