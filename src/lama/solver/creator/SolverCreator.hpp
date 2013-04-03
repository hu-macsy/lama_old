/**
 * @file SolverCreator.hpp
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
 * @brief SolverCreator.hpp
 * @author Kai Buschulte
 * @date 06.06.2012
 * $Id$
 */
#ifndef LAMA_SolverCreator_HPP_
#define LAMA_SolverCreator_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/solver/Solver.hpp>

// macros
#include <lama/macros/unique_name.hpp>

#include <lama/solver/creator/Creator.hpp>

// spirit
#include <boost/spirit/include/qi.hpp>

#include <string>

#define LAMA_PARAMETER_ENTRY( key, grammarInput, instance, className, setterFunction )             \
    -( lit(key)                                                                                       \
       >> lit('=')                                                                                       \
       >> grammarInput                                                                                   \
       [ phoenix::bind(&className::setterFunction, *dynamic_cast_<className*>( instance ), _1 ) ] \
       >> lit(';') )

namespace lama
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

class MetaSolver;

/**
 * SolverCreator: base of the Solver Registry
 */
class LAMA_DLL_IMPORTEXPORT SolverCreator: Creator
{
public:
    typedef qi::rule<std::string::const_iterator,Solver*(),ascii::space_type> RuleType;

    virtual ~SolverCreator();

    /**
     * @brief Returns the main creator rule of this solver
     */
    virtual RuleType& getCreatorRule() = 0;

protected:
    typedef qi::rule<std::string::const_iterator,void( Solver* ),ascii::space_type> InternRuleType;

    using Creator::mRScalar;
    using Creator::mRId;

    SolverCreator();

    /**
     * @brief Rule that returns SolverPtr to a linked Solver
     * this is used for preconditioning or smoother/coarse-level-solver-definitions
     */
    qi::rule<std::string::const_iterator,SolverPtr(),ascii::space_type> mRSolverReference;

    /**
     * @brief Main Rule of each Solver. Handles Loggerdefinitions.
     */
    InternRuleType mRSolver;

private:
    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief Registrates the given type.
 */
#define LAMA_SOLVERCREATOR_REGISTRATION(creatorType)                                     \
    static lama::creatorType                           \
    LAMA_UNIQUE_NAME( solRegObj, creatorType )(creatorType::id())

}

#endif // LAMA_SolverCreator_HPP_
