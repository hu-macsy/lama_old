/**
 * @file CriteriaCreator.hpp
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
 * @brief CriteriaCreator.hpp
 * @author Kai Buschulte
 * @date 29.06.2012
 * $Id$
 */

#ifndef LAMA_CRITERIACREATOR_HPP_
#define LAMA_CRITERIACREATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>

#include <lama/solver/creator/Creator.hpp>

//Norm
#include <lama/norm/L1Norm.hpp>
#include <lama/norm/L2Norm.hpp>
#include <lama/norm/MaxNorm.hpp>

//Stopping Criteria
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/criteria/ResidualStagnation.hpp>

// logging
#include <logging/logging.hpp>

#include <ostream>

namespace lama
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

class LAMA_DLL_IMPORTEXPORT CriteriaCreator: Creator
{
public:
    typedef qi::rule<std::string::const_iterator,CriterionPtr(),ascii::space_type> RuleType;
    typedef qi::symbols<char,CriterionPtr> CriteriaInstanceMap;

    /**
     * @brief Returns the lazy created singleton instance of the boolean condition creator
     */
    static CriteriaCreator& getInstance();

    /*
     * @brief Returns the grammar rule for "in place" criteria creation
     * This is used by the IterativeSolverCreator class
     * In case of i.e. criteria = IterationCount(10) AND Resi...
     * the parameter definition of a criteria is solver bound and not reusable (registered) for
     * other solver instances.
     *
     * @return  A rule that returns a BooleanContionPtr to the root of the criteria tree
     */
    static RuleType& getSolverBoundRule();

    /*
     * @brief Returns the grammar rule for independent and registered criteria creation
     * This is used by the MetaSolver class
     * In case of i.e. Criterion criteria1 = IterationCount(10) AND Resi...
     * the parameter definition of a criteria is independent and reusable (registered) so that
     * other solver instances can use this instance multiple times.
     *
     * @return  A rule that returns a BooleanContionPtr to the root of the criteria tree
     */
    static qi::rule<std::string::const_iterator,void(),ascii::space_type>& getIndependentRule();

    /**
     * @brief Returns the criteria instance map that solver rules can bind these instances
     *
     * @return Returns the criteria instance map
     */
    const CriteriaInstanceMap& getCriteriaInstanceMap();

protected:
    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:
    CriteriaCreator();

    /**
     * @brief Add/register a stopping criteria instance with a given name/id to the Factory
     *
     * @param[in] name      The name to identify the instance
     * @param[in] criteria  The instance that will be registered
     */
    void addCriteria( const std::string& name, CriterionPtr criteria );

    RuleType mRLeaf;
    RuleType mRNode;
    RuleType mRSolverBoundCriteria;
    qi::rule<std::string::const_iterator, void(), ascii::space_type> mRIndependentCriteria;
    qi::rule<std::string::const_iterator, IterationCount*(), ascii::space_type> mRIterationCount;
    qi::rule<std::string::const_iterator, ResidualThreshold*(), ascii::space_type> mRResidualThreshold;
    qi::rule<std::string::const_iterator, ResidualStagnation*(), ascii::space_type> mRResidualStagnation;
    qi::rule<std::string::const_iterator, NormPtr(), ascii::space_type> mRNormType;

    CriteriaInstanceMap mCriteriaInstanceMap;

    qi::symbols<char, ResidualThreshold::ResidualThresholdCheckMode> mResidualCheckMode;
    qi::symbols<char, Criterion::BooleanOperator> mLogicalConnective;
};

} // namespace lama

#endif // LAMA_CRITERIACREATOR_HPP_
