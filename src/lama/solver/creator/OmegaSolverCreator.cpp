/**
 * @file OmegaSolverCreator.cpp
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
 * @brief OmegaSolverCreator.cpp
 * @author Kai Buschulte
 * @date 11.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/OmegaSolverCreator.hpp>

// spirit
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/dynamic_cast.hpp>
#include <boost/spirit/home/phoenix/operator/member.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>

#include<lama/solver/DefaultJacobi.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OmegaSolverCreator::logger, "Solver.SolverCreator.OmegaSolverCreator" );

OmegaSolverCreator::OmegaSolverCreator()
    : IterativeSolverCreator()
{
    using qi::lit;
    using qi::_1;
    using qi::_a;
    using qi::_r1;
    using qi::double_;
    using qi::_val;
    using phoenix::dynamic_cast_;
    using phoenix::construct;
    using qi::debug;

    mROmegaSolver = mRIterativeSolver( _r1 ) >> LAMA_PARAMETER_ENTRY( "omega", mRScalar, _r1, OmegaSolver, setOmega );

    mROmegaSolver.name( "OmegaSolver" );
}

OmegaSolverCreator::~OmegaSolverCreator()
{
}

} //namespace lama
