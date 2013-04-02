/**
 * @file OmegaSolverCreator.hpp
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
 * @brief OmegaSolverCreator.hpp
 * @author Kai Buschulte
 * @date 11.06.2012
 * $Id$
 */
#ifndef LAMA_OmegaSolverCreator_HPP_
#define LAMA_OmegaSolverCreator_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/solver/OmegaSolver.hpp>

// base classes
#include <lama/solver/creator/IterativeSolverCreator.hpp>

//qi, phoenix, fusion
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>

namespace lama
{

namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class LAMA_DLL_IMPORTEXPORT OmegaSolverCreator: public IterativeSolverCreator
{
public:
    using IterativeSolverCreator::RuleType;

    static const std::string& id();

    virtual ~OmegaSolverCreator();

    virtual RuleType& getCreatorRule() = 0;

    typedef void (lama::OmegaSolver::*OmegaSolverFunc)( Scalar );
protected:
    OmegaSolverCreator();

    using SolverCreator::InternRuleType;
    using SolverCreator::mRScalar;

    qi::rule<std::string::const_iterator,void( Solver* ),ascii::space_type> mROmegaSolver;

    using IterativeSolverCreator::mRIterativeSolver;

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_OmegaSolverCreator_HPP_
