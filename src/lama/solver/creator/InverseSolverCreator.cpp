/**
 * @file InverseSolverCreator.cpp
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
 * @brief InverseSolverCreator.cpp
 * @author Kai Buschulte
 * @date 20.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/InverseSolverCreator.hpp>

// others
#include <lama/solver/InverseSolver.hpp>
#include <lama/solver/creator/SolverFactory.hpp>

// boost
#include <boost/config/warning_disable.hpp>
// spirit
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( InverseSolverCreator::logger, "Solver.SolverCreator.InverseSolverCreator" );

const std::string& InverseSolverCreator::id()
{
    static const std::string id = "InverseSolver";
    return id;
}

InverseSolverCreator::RuleType& InverseSolverCreator::getCreatorRule()
{
    return mRInverseSolver;
}

InverseSolverCreator::InverseSolverCreator( const std::string type )
    : SolverCreator()
{
    using qi::_r1;
    using qi::on_error;
    using qi::debug;
    using qi::fail;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_val;
    using qi::lit;

    using phoenix::val;
    using phoenix::construct;
    using phoenix::new_;

    mRInverseSolver = mRId[_val = new_<InverseSolver>( _1 )] > lit( '{' ) > mRSolver( _val ) > lit( '}' );

    mRInverseSolver.name( "InverseSolver" );

    on_error<fail>( mRInverseSolver, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );

    SolverFactory::getFactory().addSolverCreator( type, mRInverseSolver );
}

InverseSolverCreator::~InverseSolverCreator()
{
}

LAMA_SOLVERCREATOR_REGISTRATION( InverseSolverCreator );

} //namespace lama
