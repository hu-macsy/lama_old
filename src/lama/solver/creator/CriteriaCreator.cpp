/**
 * @file CriteriaCreator.cpp
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
 * @brief CriteriaCreator.cpp
 * @author Kai Buschulte
 * @date 29.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/CriteriaCreator.hpp>

// spirit
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CriteriaCreator::logger, "Solver.SolverCreator.CriteriaCreator" );

CriteriaCreator::CriteriaCreator()
    : Creator()
{
    using qi::lit;
    using qi::_val;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_a;

    using qi::on_error;
    using qi::fail;
    using qi::debug;

    using qi::double_;
    using qi::int_;
    using qi::char_;

    using phoenix::val;
    using phoenix::construct;
    using phoenix::new_;
    using phoenix::bind;
    using phoenix::ref;

    mResidualCheckMode.add( "Relative", ResidualThreshold::Relative )( "Absolute", ResidualThreshold::Absolute );

    mLogicalConnective.add( "&&", Criterion::AND )( "AND", Criterion::AND )( "OR", Criterion::OR )( "||",
            Criterion::OR );

    mRNormType = lit( "L1Norm" )[_val = construct<NormPtr>( new_<L1Norm>() )] | lit( "L2Norm" )[_val =
                     construct<NormPtr>( new_<L2Norm>() )]
                 | lit( "MaxNorm" )[_val = construct<NormPtr>( new_<MaxNorm>() )];

    mRIterationCount = lit( "IterationCount" ) > '(' >> int_[_val = new_<IterationCount>( _1 )] >> lit( ')' );

    mRResidualThreshold = lit( "ResidualThreshold" ) > '('
                          > ( mRNormType > ',' > double_ > ',' > mResidualCheckMode )[_val = new_<ResidualThreshold>( _1, _2,
                                  _3 )]
                          > lit( ')' );

    mRResidualStagnation = lit( "ResidualStagnation" ) > '('
                           > ( mRNormType > ',' > int_ > ',' > mRScalar )[_val = new_<ResidualStagnation>( _1, _2, _3 )] > ')';

    //wraps the concrete in a CriterionPtr
    mRLeaf = mRIterationCount[_val = construct<CriterionPtr>( _1 )] | mRResidualThreshold[_val =
                 construct<CriterionPtr>( _1 )] | mRResidualStagnation[_val = construct<CriterionPtr>( _1 )];

    mRNode = ( lit( '(' )
               > ( mRNode > mLogicalConnective > mRNode )[_val = construct<CriterionPtr>(
                           new_<Criterion>( _1, _3, _2 ) )] > ')' ) | mRLeaf[_val = _1];

    mRSolverBoundCriteria = ( mCriteriaInstanceMap[_val = _1] | mRNode[_val = _1] );

    mRIndependentCriteria = "Criterion"
                            > ( mRId > '=' > mRNode )[phoenix::bind( &CriteriaCreator::addCriteria, *this, _1, _2 )]
                            > lit( ';' );

    mRSolverBoundCriteria.name( "criteria" );
    mRIndependentCriteria.name( "criterion" );
    mRNode.name( "node" );
    mRLeaf.name( "leaf" );
    mRIterationCount.name( "ItCount" );
    mRResidualThreshold.name( "ResThreshold" );
    mRResidualStagnation.name( "ResStagnation" );
    mRNormType.name( "NormType" );

    on_error<fail>( mRSolverBoundCriteria, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );

    on_error<fail>( mRIndependentCriteria, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );
}

CriteriaCreator& CriteriaCreator::getInstance()
{
    static std::auto_ptr<CriteriaCreator> instance;

    if ( !instance.get() )
    {
        instance = std::auto_ptr<CriteriaCreator>( new CriteriaCreator() );
    }

    return *instance;
}

CriteriaCreator::RuleType& CriteriaCreator::getSolverBoundRule()
{
    return getInstance().mRSolverBoundCriteria;
}

qi::rule<std::string::const_iterator,void(),ascii::space_type>& CriteriaCreator::getIndependentRule()
{
    return getInstance().mRIndependentCriteria;
}

void CriteriaCreator::addCriteria( const std::string& name, CriterionPtr criteria )
{
    LAMA_ASSERT( !mCriteriaInstanceMap.find( name ), "Criteria with id " << name << " already registered." );

    mCriteriaInstanceMap.add( name, criteria );

    LAMA_LOG_INFO( logger, "Registered criteria " << name << " describing " << *criteria );
}

const CriteriaCreator::CriteriaInstanceMap& CriteriaCreator::getCriteriaInstanceMap()
{
    return mCriteriaInstanceMap;
}

} // namespace lama
