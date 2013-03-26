/**
 * @file MetaSolver.cpp
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
 * SOFTWARE.using phoenix::val;
 * @endlicense
 *
 * @brief MetaSolver.cpp
 * @author Kai Buschulte
 * @date 07.05.2012
 * $Id$
 */

// hpp
#include <lama/solver/MetaSolver.hpp>

#include <lama/solver/IterativeSolver.hpp>
#include <lama/solver/logger/LogLevel.hpp>

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
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>

#include <lama/solver/creator/CriteriaCreator.hpp>
#include <lama/solver/creator/LoggerCreator.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace lama
{

LAMA_LOG_DEF_LOGGER( MetaSolver::logger, "Solver.MetaSolver" )

MetaSolver::MetaSolver( const std::string& id )
    : Solver( id )
{
}

MetaSolver::MetaSolver( const std::string& id, LoggerPtr logger )
    : Solver( id, logger )
{
}

MetaSolver::MetaSolver( const std::string& id, const std::string& arg )
    : Solver( id )
{
    interpreteArgument( arg );
}

MetaSolver::MetaSolver( const std::string& id, const std::string& arg, LoggerPtr logger )
    : Solver( id, logger )
{
    interpreteArgument( arg );
}

MetaSolver::MetaSolverRuntime::MetaSolverRuntime()
    : SolverRuntime()
{
}

MetaSolver::~MetaSolver()
{
}

MetaSolver::MetaSolverRuntime::~MetaSolverRuntime()
{
}

void MetaSolver::initialize( const Matrix& coefficients )
{
    LAMA_LOG_INFO( logger, "Initializing MetaSolver " )

    if ( !getRuntime().mRootSolver )
    {
        LAMA_THROWEXCEPTION( "No root solver defined for MetaSolver." )
    }

    getRuntime().mRootSolver->initialize( coefficients );

    Solver::initialize( coefficients );
}

void MetaSolver::initializePreconditioner( const Matrix& coefficients, LogLevel::LogLevel level )
{
    LAMA_LOG_INFO( logger, "Initializing Meta-Preconditioner " )

    if ( !getRuntime().mRootSolver )
    {
        LAMA_THROWEXCEPTION( "No root solver defined for MetaSolver." )
    }

    IterativeSolver* rootPtr = dynamic_cast<IterativeSolver*>( &( *( getRuntime().mRootSolver ) ) );
    if ( rootPtr )
    {
        if ( rootPtr->getPreconditioner() )
        {
            rootPtr->getPreconditioner()->initialize( coefficients );
            rootPtr->getPreconditioner()->setLogLevel( level );
        }
        else
        {
            LAMA_THROWEXCEPTION( "No preconditioner defined for MetaSolver." )
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "No iterative solver defined for MetaSolver." )
    }
}

void MetaSolver::solveImpl()
{
    if ( !getRuntime().mRootSolver )
    {
        LAMA_THROWEXCEPTION( "No root solver defined for MetaSolver." )
    }

    LAMA_LOG_INFO( logger, "Solver " << getRuntime().mRootSolver->getId() << " solves now." )
    getRuntime().mRootSolver->solve( *getRuntime().mSolution, *getRuntime().mRhs );
}

void MetaSolver::interpreteArgument( const std::string arg )
{
    std::ifstream configFile;
    configFile.open( arg.c_str() );

    if ( configFile )
    {
        LAMA_LOG_DEBUG( logger, "Argument " << arg << " is a file. Reading content now." )
        std::string configuration;

        configFile.unsetf( std::ios::skipws ); // No white space skipping!
        std::copy( std::istream_iterator<char>( configFile ), std::istream_iterator<char>(),
                   std::back_inserter( configuration ) );

        parseConfiguration( configuration );
    }
    else
    {
        parseConfiguration( arg );
    }
}

void MetaSolver::parseConfiguration( const std::string& arg )
{
    LAMA_LOG_INFO( logger, "Parsing configuration " << arg )

    mLogger->startTimer( "ConfigurationTimer" );

    SolverConfigGrammar configReader; // Our grammar

//    SkipGrammar<std::string::const_iterator> skipper;

    StringIterator first = arg.begin();
    StringIterator last = arg.end();

    bool r = phrase_parse( first, last, configReader, ascii::space );

    getRuntime().mRootSolver = configReader.getRootSolver();

    LAMA_ASSERT( getRuntime().mRootSolver, "No root solver defined in this configuration." )

    LAMA_LOG_DEBUG( logger, "Solver " << *getRuntime().mRootSolver << " is root now." )

    if ( !r || first != last )
    {
        std::string rest( first, last );
        LAMA_THROWEXCEPTION( "Parsing failure. Stopped at " << rest )
    }
    mLogger->stopTimer( "ConfigurationTimer" );
    mLogger->logTime( "ConfigurationTimer", LogLevel::completeInformation,
                      "Configuration of the MetaSolver took [s]: " );
    mLogger->stopAndResetTimer( "ConfigurationTimer" );

    mLogger->logMessage( LogLevel::solverInformation, "MetaSolver configured.\n" );
}

SolverPtr MetaSolver::copy()
{
    LAMA_THROWEXCEPTION( "Copy of MetaSolver not implemented." )
    return SolverPtr();
}

MetaSolver::MetaSolverRuntime& MetaSolver::getRuntime()
{
    return mMetaSovlerRuntime;
}

const MetaSolver::MetaSolverRuntime& MetaSolver::getConstRuntime() const
{
    return mMetaSovlerRuntime;
}

LAMA_LOG_DEF_LOGGER( SolverConfigGrammar::logger, "Solver.MetaSolver.Grammar" )

SolverConfigGrammar::SolverConfigGrammar()
    : base_type( mRConfigurationSequence )
{
    using qi::lit;
    using qi::_a;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_r1;
    using qi::_val;
    using qi::lexeme;
    using qi::double_;
    using qi::int_;
    using qi::on_error;
    using qi::fail;
    using qi::debug;

    using ascii::char_;

    using phoenix::push_back;
    using phoenix::val;
    using phoenix::construct;
    using phoenix::if_;

    SolverFactory& factory = SolverFactory::getFactory();
    LAMA_LOG_DEBUG( logger, "Get creator rule symbol map." )
    const SolverFactory::TypeToCreatorMap& creatorMap = factory.getCreatorRuleSymbols();

    mRSolverConfiguration = creatorMap[_a = _1]
                            >> lazy( _a )[if_( phoenix::bind( &Solver::getId, _1 ) != "root" )[phoenix::bind(
                                    &SolverFactory::addSolver, phoenix::ref( factory ), _1 )].else_[phoenix::bind(
                                            &SolverConfigGrammar::setRootSolver, *this, construct<SolverPtr>( _1 ) )]];

    mRCriteriaConfiguration = CriteriaCreator::getIndependentRule();

    mRLoggerConfiguration = LoggerCreator::getIndependentRule();

    mRConfiguration = mRSolverConfiguration | mRCriteriaConfiguration | mRLoggerConfiguration;

    mRConfigurationSequence = *mRConfiguration;

    mRConfiguration.name( "Configuration" );
    mRConfigurationSequence.name( "ConfigurationSequence" );

    on_error<fail>( mRConfigurationSequence, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << phoenix::construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << phoenix::construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );
}

void SolverConfigGrammar::setRootSolver( SolverPtr solver )
{
    LAMA_ASSERT( solver, "Defined root solver is Null." )

    LAMA_LOG_DEBUG( logger, "Defined root solver " << *solver << "." )

    mRootSolver = solver;
}

SolverPtr SolverConfigGrammar::getRootSolver()
{
    LAMA_ASSERT( mRootSolver, "No root solver defined." )

    return mRootSolver;
}

//template <typename Iterator>
//SkipGrammar<Iterator>::SkipGrammar() : SkipGrammar::base_type(start)
//{
//    using qi::int_;
//    using qi::lit;
//    using qi::double_;
//    using qi::lexeme;
//    using ascii::char_;
//
//    comment = ascii::space | lexeme["/*" >> *char_ >> "*/"];
//
//    start %= comment;
//}

}//namespace lama
