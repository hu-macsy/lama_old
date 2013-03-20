/**
 * @file SolverCreator.cpp
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
 * @brief SolverCreator.cpp
 * @author kbuschulte
 * @date 08.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/SolverFactory.hpp>

// others
#include <lama/solver/creator/SolverCreator.hpp>
#include <lama/solver/creator/LoggerCreator.hpp>

// spirit
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/phoenix/object/dynamic_cast.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( SolverCreator::logger, "Solver.SolverCreator" );

SolverCreator::SolverCreator()
    : Creator()
{
    using qi::lexeme;
    using qi::_val;
    using qi::_1;
    using qi::_r1;
    using qi::lit;

    using ascii::char_;

    using phoenix::dynamic_cast_;
    using phoenix::ref;

    SolverFactory& factory = SolverFactory::getFactory();

    mRSolverReference = mRId[_val = phoenix::bind( &SolverFactory::getSolver, ref( factory ), _1 )];

    mRSolver = LAMA_PARAMETER_ENTRY( "logger", LoggerCreator::getSolverBoundRule(), _r1, Solver, setLogger );

    mRSolver.name( "Solver" );
    mRSolverReference.name( "SolverReference" );
    mRId.name( "VariableName" );

}

SolverCreator::~SolverCreator()
{
}

}
