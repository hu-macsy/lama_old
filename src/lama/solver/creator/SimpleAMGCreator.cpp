/**
 * @file SimpleAMGCreator.cpp
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
 * @brief SimpleAMGCreator.cpp
 * @author Kai Buschulte
 * @date 11.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/SimpleAMGCreator.hpp>

// others
#include <lama/solver/SimpleAMG.hpp>
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
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/spirit/home/phoenix/object/dynamic_cast.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>

#include <boost/get_pointer.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( SimpleAMGCreator::logger, "Solver.SolverCreator.SimpleAMGCreator" );

SimpleAMGCreator::SimpleAMGCreator( const std::string type )
    : IterativeSolverCreator()
{
    using qi::_r1;
    using qi::_r2;
    using qi::lit;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_a;
    using qi::int_;
    using qi::_val;

    using phoenix::ref;
    using phoenix::new_;
    using phoenix::construct;
    using phoenix::dynamic_cast_;
    using phoenix::val;

    using qi::on_error;
    using qi::fail;
    using qi::debug;

    mRSimpleAMG =
        mRId[_val = new_<SimpleAMG>( _1 )] > lit( '{' ) > mRIterativeSolver( _val )
        > LAMA_PARAMETER_ENTRY( "smoother", mRSolverReference, _val, SimpleAMG, setSmoother )
        > LAMA_PARAMETER_ENTRY( "coarseLevelSolver", mRSolverReference, _val, SimpleAMG, setCoarseLevelSolver )
        > LAMA_PARAMETER_ENTRY( "maxLevels", int_, _val, SimpleAMG, setMaxLevels )
        > LAMA_PARAMETER_ENTRY( "minVarsCoarseLevel", int_, _val, SimpleAMG, setMinVarsCoarseLevel )
        > LAMA_PARAMETER_ENTRY( "hostOnlyLevel", int_, _val, SimpleAMG, setHostOnlyLevel )
        > LAMA_PARAMETER_ENTRY( "replicatedLevel", int_, _val, SimpleAMG, setReplicatedLevel )
        > lit( '}' );

    mRSimpleAMG.name( "SimpleAMG" );

    on_error<fail>( mRSimpleAMG, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );

    SolverFactory::getFactory().addSolverCreator( type, mRSimpleAMG );
}

SimpleAMGCreator::~SimpleAMGCreator()
{
}

const std::string& SimpleAMGCreator::id()
{
    static const std::string id = "SimpleAMG";
    return id;
}

SimpleAMGCreator::RuleType& SimpleAMGCreator::getCreatorRule()
{
    return mRSimpleAMG;
}

LAMA_SOLVERCREATOR_REGISTRATION( SimpleAMGCreator );

} //namespace lama
