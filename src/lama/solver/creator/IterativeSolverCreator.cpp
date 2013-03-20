/**
 * @file IterativeSolverCreator.cpp
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
 * @brief IterativeSolverCreator.cpp
 * @author Kai Buschulte
 * @date 11.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/IterativeSolverCreator.hpp>

// others
#include <lama/solver/IterativeSolver.hpp>
#include <lama/solver/creator/CriteriaCreator.hpp>

// spirit
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/dynamic_cast.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( IterativeSolverCreator::logger, "Solver.SolverCreator.IterativeSolverCreator" );

IterativeSolverCreator::IterativeSolverCreator()
    : SolverCreator()
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    using qi::lit;
    using qi::_val;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_r1;
    using qi::_r2;
    using qi::double_;
    using qi::on_error;
    using qi::fail;
    using qi::debug;

    using phoenix::ref;
    using phoenix::construct;
    using phoenix::val;
    using phoenix::dynamic_cast_;
    using phoenix::new_;
    using phoenix::bind;

    mRIterativeSolver =
        mRSolver( _r1 ) >> LAMA_PARAMETER_ENTRY( "preconditioner", mRSolverReference, _r1, IterativeSolver, setPreconditioner)
        >> LAMA_PARAMETER_ENTRY( "stoppingCriterion", CriteriaCreator::getSolverBoundRule(), _r1, IterativeSolver, setStoppingCriterion);

    mRPreconditioner.name( "Preconditioner" );
    mRIterativeSolver.name( "IterativeSolver" );
}

IterativeSolverCreator::~IterativeSolverCreator()
{
}

} //namespace lama
